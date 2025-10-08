#!/bin/bash

# Set the Python environment (if needed)
# module load Miniconda3
conda activate htmd
module load VMD
module load Chimera

# Define paths to Python scripts
SCRIPT_DIR="pipeline"
CONFIG_SCRIPT="$SCRIPT_DIR/config_pipeline.py"
LIGAND_PARAM_SCRIPT="$SCRIPT_DIR/run_ligand_parametrization.py"
BUILD_MODEL_SCRIPT="$SCRIPT_DIR/run_build_model.py"
EQUILIBRATION_SCRIPT="$SCRIPT_DIR/run_equilibration.py"
PRODUCTION_SCRIPT="$SCRIPT_DIR/run_production.py"
WRAP_SCRIPT="$SCRIPT_DIR/run_wrap.py"
QUALITY_SCRIPT="$SCRIPT_DIR/run_quality.py"

# Function to run a step
run_step() {
    local step_number=$1
    local step_name=$2
    local script=$3

    echo ""
    echo -e "\e[34mStep $step_number: Preparing $step_name...\e[0m"
    echo ""
    local start_time=$(date +%s)
    /home/agarcia/miniconda3/envs/htmd/bin/python $script
    if [ $? -ne 0 ]; then
        echo "Error during $step_name. Exiting."
        exit 1
    fi
    local end_time=$(date +%s)
    local elapsed_time=$((end_time - start_time))
    echo ""
    local days=$((elapsed_time / 86400))
    local hours=$(( (elapsed_time % 86400) / 3600 ))
    local minutes=$(( (elapsed_time % 3600) / 60 ))
    local seconds=$((elapsed_time % 60))
    echo "Step $step_number: $step_name prepared in ${days}d ${hours}h ${minutes}m ${seconds}s."
    echo ""
}

# Helper function to display usage information
show_help() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --no-prompt          Run the pipeline without prompting for steps."
    echo "  --steps <steps>      Specify steps to run (e.g., '1 2 3')."
    echo "  --help               Show this help message and exit."
    echo ""
    echo "Steps:"
    echo "  1: Ligand Parametrization"
    echo "  2: Build Model"
    echo "  3: Equilibration"
    echo "  4: Production"
    echo "  5: Wrap"
    echo "  6: Quality Check"
    echo ""
}

echo "###################################################"
echo "Molecular Dynamics Simulation Pipeline using ACEMD"
echo "###################################################"
echo " "

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --no-prompt)
            no_prompt=true
            shift
            ;;
        --steps)
            steps=$2
            shift 2
            ;;
        --help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Step 0: Configuration (always executed)
if [ -z "$no_prompt" ]; then
    echo "Configuration script will always be executed."
    echo ""
    echo -e "\033[93mWARNING: Ensure you have checked and updated the script config_pipeline.py and inputs as needed before running the pipeline.\033[0m"
    echo -e "\033[93mIf you do not want to see this message again, use the flag --no-prompt.\033[0m"
    echo ""
    read -p "Press Enter to continue..."
fi
run_step 0 "configuration" $CONFIG_SCRIPT

# Prompt user for steps to run if --no-prompt is not set and --steps is not provided
if [[ -z "$steps" && -z "$no_prompt" ]]; then
    echo "Select steps to run (e.g., 1 2 3)."
    echo "1: Ligand Parametrization"
    echo "2: Build Model"
    echo "3: Equilibration"
    echo "4: Production"
    echo "5: Wrap"
    echo "6: Quality Check"
    read -p "Enter step numbers: " steps
fi

# Run selected steps
if [[ ! -z "$steps" ]]; then
    for step in $steps; do
        case $step in
            1) run_step 1 "ligand parametrization" $LIGAND_PARAM_SCRIPT ;;
            2) run_step 2 "build model" $BUILD_MODEL_SCRIPT ;;
            3) run_step 3 "equilibration" $EQUILIBRATION_SCRIPT ;;
            4) run_step 4 "production" $PRODUCTION_SCRIPT ;;
            5) run_step 5 "wrap" $WRAP_SCRIPT ;;
            6) run_step 6 "quality check" $QUALITY_SCRIPT ;;
            *) echo "Invalid step: $step. Skipping." ;;
        esac
    done
else
    # Default steps to run if no --steps flag is provided
    run_step 1 "ligand parametrization" $LIGAND_PARAM_SCRIPT
    run_step 2 "build model" $BUILD_MODEL_SCRIPT
fi

echo "Pipeline completed successfully!"
