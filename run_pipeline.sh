#!/bin/bash

# Set the Python environment (if needed)
module load Miniconda3
module load VMD
module load Chimera
module load NAMD
source activate pipeline/htmd

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
    shift 3
    local args=$@

    echo ""
    echo -e "\e[34mStep $step_number: Preparing $step_name...\e[0m"
    echo ""
    local start_time=$(date +%s)
    ./pipeline/htmd/bin/python $script $args
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
    echo "  --input-json <file>  Path to the input JSON file (required)."
    echo "  --config-json <file> Path to the configuration JSON file (required)."
    echo "  --no-prompt          Run the pipeline without prompting for steps."
    echo "  --steps <steps>      Specify steps to run (e.g., '1 2 3')."
    echo "  --instructions       Show instructions and exit."
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

# Helper function to display instructions
show_instructions() {
    echo "###################################################"
    echo "Instructions for Running the Molecular Dynamics Pipeline"
    echo "###################################################"
    echo ""
    echo "0. (OPTIONAL) Use the script generate_input_json.py to create a template inputs.json file."
    echo "1. Ensure you have the required input files:"
    echo "   - inputs.json: Contains the input data for the pipeline."
    echo "   - config.json: Contains the configuration settings for the pipeline."
    echo ""
    echo "2. Use the example files (inputs_empty.json and config_empty.json) as templates."
    echo "   - Modify these files to suit your specific requirements."
    echo "     Example:"
    echo ""
    echo "   - inputs.json:"
    echo "     ["
    echo "         {"
    echo "             \"name\": \"delta_opioid_naltrindole\", // A unique name for this simulation system"  
    echo "             \"pdbfile\": \"demo/4EJ4.pdb\",  // Relative to basepath to pdb. NOT ABSOLUTE PATH"
    echo "             \"modres\": [], // List of modified residues in the PDB file (e.g \"MODIFIED_RESNAME1\", \"MODIFIED_RESNAME\")" 
    echo "             \"ligands\": [{"
    echo "                             \"resname\": \"EJ4\", // Residue name of the small molecule ligand in the PDB file"
    echo "                             \"name\": \"Naltrindole\", // A unique name for this small molecule ligand"
    echo "                             \"covalently_bound\": false,  // If the SMALMOL is covalently bound to peptide"
    echo "                             \"inchikey\": \"WIYUZYBFCWCCQJ-IFKAHUTRSA-N\""
    echo "                         }],"
    echo "             \"apo\": False,  // Do you wish to simulate an apo-version of this structure, by removing all ligands and not-main-proteins of it?"
    echo "             \"prot_chain\": \"A\",  // Chain ID of the main protein in this system (a GPCR, for us)"
    echo "             \"pdbcode\": \"4EJ4\",  // Closest-ressembling PDB structure of this GPCR. Put False if there is none"
    echo "             \"curated\": false,  // If system has been already properly protonated"
    echo "             \"sod2x50\": true,  // If system requires addition of sodium near 2x50"
    echo "             \"isgpcr\": true  // If system is indeed a GPCR"
    echo "         }"
    echo "     ]"
    echo ""
    echo "   - config.json:"
    echo "   - IMPORTANT: Only include the parameters that you need, remove the not used. The rest will be set to default values."
    echo "     Example:"
    echo "     {"
    echo "         \"paramchem_credentials\": {"
    echo "             \"username\": null, // CCgenF username (required)"
    echo "             \"password\": null // CCgenF password (required)"
    echo "         },"
    echo "         \"machine_settings\": {"
    echo "             \"device_gpu\": 4, // GPU device ID to be used for the simulations"
    echo "             \"strucpath\": \"./input_structures/\", // Path to the input structures folder"
    echo "             \"resultspath\": \"./simulation_output/\" // Path to the output results folder"
    echo "         },"
    echo "         \"simulation_parameters\": {"
    echo "             \"membrane_path_pdb\": \"membrane/popc36_box_renumbered.pdb\", // Path to the membrane PDB file"
    echo "             \"topparpath\": \"toppar/TOP_PARAMS_ACE3/\", // Path to the Topology and Parameter files"
    echo "             \"toposfilenames\": ["
    echo "                 \"General_top_params/topologies/top_all36_cgenff.rtf\","
    echo "                 \"General_top_params/topologies/top_all36_prot.rtf\","
    echo "                 \"General_top_params/topologies/top_all36_na.rtf\","
    echo "                 \"General_top_params/topologies/top_all36_lipid.rtf\","
    echo "                 \"General_top_params/topologies/top_all36_carb.rtf\","
    echo "                 \"Specific_top_params/topologies/toppar_all36_lipid_ether.rtf\","
    echo "                 \"Specific_top_params/topologies/toppar_water_ions.rtf\","
    echo "                 \"Specific_top_params/topologies/toppar_all36_na_nad_ppi.rtf\","
    echo "                 \"Specific_top_params/topologies/toppar_all36_prot_na_combined.rtf\","
    echo "                 \"David_top_params/topologies/my_patches.rtf\","
    echo "                 \"David_top_params/topologies/toppar_all36_lipid_cholesterol_CLR.rtf\","
    echo "                 \"David_top_params/topologies/toppar_all36_prot_retinol.rtf\""
    echo "             ], // Topology files (.rtf)"
    echo "             \"paramsfilesnames\": ["
    echo "                 \"General_top_params/parameters/par_all36_cgenff.prm\","
    echo "                 \"General_top_params/parameters/par_all36m_prot.prm\","
    echo "                 \"General_top_params/parameters/par_all36_na.prm\","
    echo "                 \"General_top_params/parameters/par_all36_lipid.prm\","
    echo "                 \"General_top_params/parameters/par_all36_carb.prm\","
    echo "                 \"Specific_top_params/parameters/toppar_all36_lipid_cholesterol.prm\","
    echo "                 \"Specific_top_params/parameters/toppar_all36_lipid_ether.prm\","
    echo "                 \"Specific_top_params/parameters/toppar_water_ions.prm\","
    echo "                 \"Specific_top_params/parameters/toppar_all36_na_nad_ppi.prm\","
    echo "                 \"Specific_top_params/parameters/toppar_all36_prot_na_combined.prm\","
    echo "                 \"Specific_top_params/parameters/toppar_all36_prot_retinol.prm\","
    echo "                 \"David_top_params/parameters/modres_params.prm\","
    echo "                 \"David_top_params/parameters/modres_crossterm.prm\""
    echo "             ], // Parameter files (.prm)"
    echo "             \"streamfilenames\": [], // Stream files (topology+parameters)"
    echo "             \"ligandsdict_path\": \"./ligands.json\", // Path to the ligands dictionary"
    echo "             \"modres_path\": \"./modified_residues.json\", // Path to the modified residues dictionary"
    echo "             \"preparation\": {"
    echo "                 \"new_pdb_chain\": \"R\", // New chain ID to assign to the protein"
    echo "                 \"membrane_lipid_segid\": \"MEM\", // SegID to assign to membrane lipids"
    echo "                 \"coldist\": 1.3, // Distance below which two atoms are considered to collide"
    echo "                 \"membrane_distance\": 20, // Distance between the pbc box and the space to be filled with membrane atoms"
    echo "                 \"water_thickness\": 20, // Size in Z-axis of the solvation water layers"
    echo "                 \"buffer\": 2.4, // Distance between solvation waters and the rest of the system"
    echo "                 \"water_margin\": 4, // Distance in the Z-axis to be penetrated by the solvation box to avoid the formation of a VOID"
    echo "                 \"lipidlike_blacklist\": [\"OLA\", \"OLB\", \"OLC\", \"PLM\", \"HTG\", \"LPP\", \"PEF\", \"2CV\", \"SOG\", \"TWT\", \"STE\", \"LMT\", \"MPG\", \"DGA\", \"1WV\", \"POV\", \"FLC\", \"PGW\", \"PC1\", \"LDA\", \"J40\", \"BNG\"], // Lipid-like residues to ignore"
    echo "                 \"detergent_blacklist\": [\"OLA\", \"OLB\", \"PEG\", \"GOL\", \"BOG\", \"OLC\", \"P6G\", \"P33\", \"UNL\", \"UNX\", \"PLM\", \"HTG\", \"12P\", \"LPP\", \"PEF\", \"2CV\", \"SOG\", \"TWT\", \"PGE\", \"SO4\", \"STE\", \"LMT\", \"ACT\", \"ACE\", \"MHA\", \"CIT\", \"1PE\", \"MPG\", \"EPE\", \"PG4\", \"DGA\", \"PO4\", \"DMS\", \"TAR\", \"1WV\", \"EDO\", \"BU1\", \"ACM\", \"PG6\", \"TLA\", \"SCN\", \"TCE\", \"MES\", \"EDT\", \"POV\", \"MLI\", \"SIN\", \"PGO\", \"FLC\", \"HTO\", \"PGW\", \"NO3\", \"PE5\", \"NH2\", \"NME\", \"NH4\", \"IOD\", \"ZN\", \"BR\", \"HG\", \"PC1\", \"LDA\", \"C8E\", \"J40\", \"BNG\"], // Detergent residues to ignore"
    echo "                 \"glucids_blacklist\": [\"MAN\", \"NAG\", \"BGC\", \"TRE\", \"9NS\", \"BMA\", \"FUC\", \"A2G\"], // Carbohydrate residues to ignore"
    echo "                 \"noparams_ligs\": [\"RET\", \"CLR\", \"LYR\", \"ALA\", \"ARG\", \"ASN\", \"ASP\", \"CYS\", \"GLU\", \"GLN\", \"GLY\", \"HIS\", \"ILE\", \"LEU\", \"LYS\", \"MET\", \"PHE\", \"PRO\", \"SER\", \"THR\", \"TRP\", \"TYR\", \"VAL\", \"AMP\", \"ADP\", \"ATP\", \"GMP\", \"GDP\", \"GTP\", \"TMP\", \"TDP\", \"TTP\", \"CMP\", \"CDP\", \"CTP\", \"SEP\", \"TPO\"] // Ligands for which parameters are not required"
    echo "             },"
    echo "             \"equilibration\": {"
    echo "                 \"temperature\": 310, // Target temperature for equilibration"
    echo "                 \"minim_steps\": 5000, // Number of minimization steps before equilibration"
    echo "                 \"equil_timestep\": 2, // Timestep to be used during equilibration"
    echo "                 \"equil_simtime\": 40, // Total simulation time for equilibration"
    echo "                 \"const_sel\": \"protein and name C CA N O or not (protein or lipid or water or ions) and noh or segid IO WAT and noh\" // Selection of atoms to be constrained during equilibration"
    echo "             },"
    echo "             \"production\": {"
    echo "                 \"timestep\": 4, // Simulation timestep to be used during production (femtoseconds)"
    echo "                 \"trajperiod\": 50000, // Simulation period: number of timesteps occurring between frames"
    echo "                 \"repnum\": 3, // Number of replicas to be run during production"
    echo "                 \"prod_simtime\": 500 // Total simulation time for production"
    echo "             }"
    echo "         }"
    echo "     }"
    echo ""
    echo "3. Run the pipeline using the following command:"
    echo "   ./run_pipeline.sh --input-json <path_to_inputs.json> --config-json <path_to_config.json>"
    echo ""
    echo "4. Use additional flags as needed:"
    echo "   - --no-prompt: Run the pipeline without interactive prompts."
    echo "   - --steps: Specify the steps to run (e.g., '1 2 3')."
    echo ""
    echo "Steps:"
    echo "1: Ligand Parametrization"
    echo "   - This step involves parametrizing the ligand molecules for the simulation."
    echo "2: Build Model"
    echo "   - This step builds the molecular model, including proteins, ligands, and membranes."
    echo "3: Equilibration"
    echo "   - This step equilibrates the system to prepare it for production runs."
    echo "4: Production"
    echo "   - This step performs the production molecular dynamics simulations."
    echo "5: Wrap"
    echo "   - This step wraps the trajectory files for visualization and analysis."
    echo "6: Quality Check"
    echo "   - This step checks the quality of the simulation results."
    echo ""
    echo "5. For more information, use the --help flag."
    echo ""
    echo "###################################################"
    exit 0
}

echo "###################################################"
echo "Molecular Dynamics Simulation Pipeline using ACEMD"
echo "###################################################"
echo " "

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input-json)
            input_json=$2
            shift 2
            ;;
        --config-json)
            config_json=$2
            shift 2
            ;;
        --no-prompt)
            no_prompt=true
            shift
            ;;
        --steps)
            steps=$2
            shift 2
            ;;
        --instructions)
            show_instructions
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

# Ensure required parameters are provided
if [[ -z "$input_json" || -z "$config_json" ]]; then
    echo "Error: --input-json and --config-json are required."
    show_help
    exit 1
fi

# Step 0: Configuration (always executed)
if [ -z "$no_prompt" ]; then
    echo "Configuration script will always be executed."
    echo ""
    echo -e "\033[93mWARNING: Ensure you have checked and updated the inputs (inputs.json & config.json) as needed before running the pipeline.\033[0m"
    echo -e "\033[93mYou have example files in the pipeline directory: inputs_empty.json & config_empty.json. Use the flag --instructions to see both examples and explanation how to run this pipeline.\033[0m"
    echo -e "\033[93mIf you do not want to see this message again, use the flag --no-prompt.\033[0m"
    echo ""
    read -p "Press Enter to continue..."
fi
run_step 0 "configuration" $CONFIG_SCRIPT --input-json $input_json --config-json $config_json

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
            1) run_step 1 "ligand parametrization" $LIGAND_PARAM_SCRIPT --input-json $input_json --config-json $config_json;;
            2) run_step 2 "build model" $BUILD_MODEL_SCRIPT --input-json $input_json --config-json $config_json;;
            3) run_step 3 "equilibration" $EQUILIBRATION_SCRIPT --input-json $input_json --config-json $config_json;;
            4) run_step 4 "production" $PRODUCTION_SCRIPT --input-json $input_json --config-json $config_json;;
            5) run_step 5 "wrap" $WRAP_SCRIPT --input-json $input_json --config-json $config_json;;
            6) run_step 6 "quality check" $QUALITY_SCRIPT --input-json $input_json --config-json $config_json;;
            *) echo "Invalid step: $step. Skipping." ;;
        esac
    done
else
    # Default steps to run if no --steps flag is provided
    run_step 1 "ligand parametrization" $LIGAND_PARAM_SCRIPT --input-json $input_json --config-json $config_json
    run_step 2 "build model" $BUILD_MODEL_SCRIPT --input-json $input_json --config-json $config_json
fi

echo "Pipeline completed successfully!"
