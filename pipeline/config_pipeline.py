"""
Simulation Pipeline Configuration Script

This script is part of a simulation pipeline designed to prepare and execute molecular dynamics simulations.
It provides configuration settings, paths, and user inputs required for the pipeline to function properly.

Key Features:
- Defines paths to required software, input files, and output directories.
- Loads topology, parameter, and stream files for molecular simulations.
- Allows user customization through a JSON dictionary for system-specific data.
- Supports integration with external tools like Chimera, VMD, and SLURM.

Ensure all paths and credentials are correctly set before running the pipeline.

# USER INPUTS: introduce two JSON dictionary, one with the data for your system in the following format (input.json):

Example:
[
    {
        "name": "NAME",
        "pdbfile": "PATH/TO/YOUR_PDB.pdb",  # Relative to basepath. NOT ABSOLUTE PATH
        "modres": ["MODIFIED_RESNAME1", "MODIFIED_RESNAME"],
        "ligands": [{
                        "resname": "SMALMOL_RESNAME1",
                        "name": "SMALMOL_NAME1",
                        "covalently_bound": true,  # If the SMALMOL is covalently bound to peptide
                        "inchikey": "RLDFVSFNVOLFEU-UHFFFAOYSA-N"
                    }],
        "apo": False,  # Do you wish to simulate an apo-version of this structure, by removing all ligands and not-main-proteins of it?
        "prot_chain": "R",  # Chain ID of the main protein in this system (a GPCR, for us)
        "pdbcode": "PDBCODE",  # Closest-ressembling PDB structure of this GPCR. Put False if there is none
        "curated": true,  # If system has been already properly protonated
        "sod2x50": true,  # If system requires addition of sodium near 2x50
        "isgpcr": true  # If system is indeed a GPCR
    }
]

# Another JSON dictionary with user configuration for the pipeline (config.json):

{
    "paramchem_credentials": {
        "username": null, // CCgenF username 
        "password": null // CCgenF password
    },
    "machine_settings": {
        "device_gpu": 4, // GPU device ID to be used for the simulations
        "strucpath": "./input_structures/", // Path to the input structures folder
        "resultspath": "./simulation_output/" // Path to the output results folder
    },
    "simulation_parameters": { 
        "membrane_path_pdb": "membrane/popc36_box_renumbered.pdb", // Path to the membrane PDB file
        "topparpath": "toppar/TOP_PARAMS_ACE3/", // Path to the Topology and Parameter files
        "toposfilenames": [
            "General_top_params/topologies/top_all36_cgenff.rtf",
            "General_top_params/topologies/top_all36_prot.rtf",
            "General_top_params/topologies/top_all36_na.rtf",
            "General_top_params/topologies/top_all36_lipid.rtf",
            "General_top_params/topologies/top_all36_carb.rtf",
            "Specific_top_params/topologies/toppar_all36_lipid_ether.rtf",
            "Specific_top_params/topologies/toppar_water_ions.rtf",
            "Specific_top_params/topologies/toppar_all36_na_nad_ppi.rtf",
            "Specific_top_params/topologies/toppar_all36_prot_na_combined.rtf",
            "David_top_params/topologies/my_patches.rtf",
            "David_top_params/topologies/toppar_all36_lipid_cholesterol_CLR.rtf",
            "David_top_params/topologies/toppar_all36_prot_retinol.rtf"
        ], // Topology files (.rtf)
        "paramsfilesnames": [
            "General_top_params/parameters/par_all36_cgenff.prm",
            "General_top_params/parameters/par_all36m_prot.prm",
            "General_top_params/parameters/par_all36_na.prm",
            "General_top_params/parameters/par_all36_lipid.prm",
            "General_top_params/parameters/par_all36_carb.prm",
            "Specific_top_params/parameters/toppar_all36_lipid_cholesterol.prm",
            "Specific_top_params/parameters/toppar_all36_lipid_ether.prm",
            "Specific_top_params/parameters/toppar_water_ions.prm",
            "Specific_top_params/parameters/toppar_all36_na_nad_ppi.prm",
            "Specific_top_params/parameters/toppar_all36_prot_na_combined.prm",
            "Specific_top_params/parameters/toppar_all36_prot_retinol.prm",
            "David_top_params/parameters/modres_params.prm",
            "David_top_params/parameters/modres_crossterm.prm"
        ], // Parameter files (.prm)
        "streamfilenames": [], // Stream files (topology+parameters)
        "ligandsdict_path": "./ligands.json", // Path to the ligands dictionary
        "modresdict_path": "./modified_residues.json", // Path to the modified residues dictionary
        "preparation": {
            "new_pdb_chain": "R", // New chain ID to assign to the protein
            "membrane_lipid_segid": "MEM", // SegID to assign to membrane lipids
            "ph": 7.4, // PH value for protonation
            "remove_sod2x50": False, // remove the sodium 2x50 and let the residue protonated
            "coldist": 1.3, // Distance below which two atoms are considered to collide
            "membrane_distance": 20, // Distance between the pbc box and the space to be filled with membrane atoms
            "water_thickness": 20, // Size in Z-axis of the solvation water layers
            "buffer": 2.4, // Distance between solvation waters and the rest of the system
            "water_margin": 4, // Distance in the Z-axis to be penetrated by the solvation box to avoid the formation of a VOID
            "lipidlike_blacklist": ["OLA", "OLB", "OLC", "PLM", "HTG", "LPP", "PEF", "2CV", "SOG", "TWT", "STE", "LMT", "MPG", "DGA", "1WV", "POV", "FLC", "PGW", "PC1", "LDA", "J40", "BNG"], // Lipid-like residues to ignore
            "detergent_blacklist": ["OLA", "OLB", "PEG", "GOL", "BOG", "OLC", "P6G", "P33", "UNL", "UNX", "PLM", "HTG", "12P", "LPP", "PEF", "2CV", "SOG", "TWT", "PGE", "SO4", "STE", "LMT", "ACT", "ACE", "MHA", "CIT", "1PE", "MPG", "EPE", "PG4", "DGA", "PO4", "DMS", "TAR", "1WV", "EDO", "BU1", "ACM", "PG6", "TLA", "SCN", "TCE", "MES", "EDT", "POV", "MLI", "SIN", "PGO", "FLC", "HTO", "PGW", "NO3", "PE5", "NH2", "NME", "NH4", "IOD", "ZN", "BR", "HG", "PC1", "LDA", "C8E", "J40", "BNG"], // Detergent residues to ignore
            "glucids_blacklist": ["MAN", "NAG", "BGC", "TRE", "9NS", "BMA", "FUC", "A2G"], // Carbohydrate residues to ignore
            "noparams_ligs": ["RET", "CLR", "LYR", "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "AMP", "ADP", "ATP", "GMP", "GDP", "GTP", "TMP", "TDP", "TTP", "CMP", "CDP", "CTP", "SEP", "TPO"], // Ligands for which parameters are not required
            "tym_remove": False // Whether to remove TYMs from the structure
        },
        "equilibration": {
            "temperature": 310, // Target temperature for equilibration
            "minim_steps": 5000, // Number of minimization steps before equilibration
            "equil_timestep": 2, // Timestep to be used during equilibration
            "equil_simtime": 40, // Total simulation time for equilibration
            "const_sel": "protein and name C CA N O or not (protein or lipid or water or ions) and noh or segid IO WAT and noh" // Selection of atoms to be constrained during equilibration
        },
        "production": {
            "timestep": 4, // Simulation timestep to be used during production (femtoseconds)
            "trajperiod": 50000, // Simulation period: number of timesteps occurring between frames
            "repnum": 3, // Number of replicas to be run during production
            "prod_simtime": 500 // Total simulation time for production
        }
    }
}
"""

import sys
import json
import argparse
import os
import subprocess

# Define default values
DEFAULTS = {
        "paramchem_credentials": {
            "parametrization": True,
            "username": "",
            "password": ""
        },
        "machine_settings": {
            "device_gpu": 4,
            "strucpath": "./input_structures/",
            "resultspath": "./simulation_output/"
        },
        "simulation_parameters": { 
            "membrane_path_pdb": "membrane/popc36_box_renumbered.pdb",
            "topparpath": "toppar/TOP_PARAMS_ACE3/",
            "toposfilenames": [
                "General_top_params/topologies/top_all36_cgenff.rtf",
                "General_top_params/topologies/top_all36_prot.rtf",
                "General_top_params/topologies/top_all36_na.rtf",
                "General_top_params/topologies/top_all36_lipid.rtf",
                "General_top_params/topologies/top_all36_carb.rtf",
                "Specific_top_params/topologies/toppar_all36_lipid_ether.rtf",
                "Specific_top_params/topologies/toppar_water_ions.rtf",
                "Specific_top_params/topologies/toppar_all36_na_nad_ppi.rtf",
                "Specific_top_params/topologies/toppar_all36_prot_na_combined.rtf",
                "David_top_params/topologies/my_patches.rtf",
                "David_top_params/topologies/toppar_all36_lipid_cholesterol_CLR.rtf",
                "David_top_params/topologies/toppar_all36_prot_retinol.rtf"
            ],
            "paramsfilesnames": [
                "General_top_params/parameters/par_all36_cgenff.prm",
                "General_top_params/parameters/par_all36m_prot.prm",
                "General_top_params/parameters/par_all36_na.prm",
                "General_top_params/parameters/par_all36_lipid.prm",
                "General_top_params/parameters/par_all36_carb.prm",
                "Specific_top_params/parameters/toppar_all36_lipid_cholesterol.prm",
                "Specific_top_params/parameters/toppar_all36_lipid_ether.prm",
                "Specific_top_params/parameters/toppar_water_ions.prm",
                "Specific_top_params/parameters/toppar_all36_na_nad_ppi.prm",
                "Specific_top_params/parameters/toppar_all36_prot_na_combined.prm",
                "Specific_top_params/parameters/toppar_all36_prot_retinol.prm",
                "David_top_params/parameters/modres_params.prm",
                "David_top_params/parameters/modres_crossterm.prm"
            ],
            "streamfilenames": [],
            "ligandsdict_path": "./ligands.json",
            "modres_path": "./modified_residues.json",
            "preparation": {
                "new_pdb_chain": "R",
                "membrane_lipid_segid": "MEM",
                "ph": 7.4,
                "remove_sod2x50": False,
                "coldist": 1.3,
                "membrane_distance": 20,
                "water_thickness": 20,
                "buffer": 2.4,
                "water_margin": 4,
                "lipidlike_blacklist": ["OLA", "OLB", "OLC", "PLM", "HTG", "LPP", "PEF", "2CV", "SOG", "TWT", "STE", "LMT", "MPG", "DGA", "1WV", "POV", "FLC", "PGW", "PC1", "LDA", "J40", "BNG"],
                "detergent_blacklist": ["OLA", "OLB", "PEG", "GOL", "BOG", "OLC", "P6G", "P33", "UNL", "UNX", "PLM", "HTG", "12P", "LPP", "PEF", "2CV", "SOG", "TWT", "PGE", "SO4", "STE", "LMT", "ACT", "ACE", "MHA", "CIT", "1PE", "MPG", "EPE", "PG4", "DGA", "PO4", "DMS", "TAR", "1WV", "EDO", "BU1", "ACM", "PG6", "TLA", "SCN", "TCE", "MES", "EDT", "POV", "MLI", "SIN", "PGO", "FLC", "HTO", "PGW", "NO3", "PE5", "NH2", "NME", "NH4", "IOD", "ZN", "BR", "HG", "PC1", "LDA", "C8E", "J40", "BNG"],
                "glucids_blacklist": ["MAN", "NAG", "BGC", "TRE", "9NS", "BMA", "FUC", "A2G"],
                "noparams_ligs": ["RET", "CLR", "LYR", "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "AMP", "ADP", "ATP", "GMP", "GDP", "GTP", "TMP", "TDP", "TTP", "CMP", "CDP", "CTP", "SEP", "TPO"],
                "tym_remove": False
            },
            "equilibration": {
                "temperature": 310,
                "minim_steps": 5000,
                "equil_timestep": 2,
                "equil_simtime": 40,
                "const_sel": "protein and name C CA N O or not (protein or lipid or water or ions) and noh or segid IO WAT and noh"
            },
            "production": {
                "timestep": 4,
                "trajperiod": 50000,
                "repnum": 3,
                "prod_simtime": 500
            }
        }
    }

# Add Chimera's and VMD paths (replace by your own)
def detect_path(executable):
    try:
        path = subprocess.check_output(['which', executable], stderr=subprocess.DEVNULL).decode().strip()
        print(f"Detected path for {executable}: {path}")
        return path
    except subprocess.CalledProcessError:
        print(f"Executable {executable} not found in PATH.")
        return None

def merge_with_defaults(input_data, defaults):
    """
    Recursively merge input data with default values.
    If a key is missing or has a null/empty value, the default value is used.
    """
    for key, default_value in defaults.items():
        if key not in input_data or input_data[key] in [None, "", [], {}]:
            input_data[key] = default_value
        elif isinstance(default_value, dict):
            # Recursively merge nested dictionaries
            input_data[key] = merge_with_defaults(input_data.get(key, {}), default_value)
    return input_data

#if __name__ == "__main__":
# Parse command-line arguments
parser = argparse.ArgumentParser(description="Update configuration JSON with default values.")
parser.add_argument("--input-json", required=True, help="Path to the input JSON file.")
parser.add_argument("--config-json", required=True, help="Path to the output JSON file.")
args = parser.parse_args()

# Load the config input JSON file
with open(args.config_json, "r") as f:
    input_data = json.load(f)

# Merge with default values
updated_data = merge_with_defaults(input_data, DEFAULTS)

# Passwords for paramchem
parametrization = updated_data['paramchem_credentials']['parametrization']
username = updated_data['paramchem_credentials']['username']
password = updated_data['paramchem_credentials']['password']
if parametrization:
    if not (username and password) or (username == "" or password == ""):
        raise Exception("Please define your password and username for CGenFF app on config_pipeline.py. It is required to parameterize small molecules")

# Paths
basepath = os.getcwd() + '/'

scriptspath = basepath + 'pipeline/'
strucpath = updated_data['machine_settings']['strucpath']
resultspath = updated_data['machine_settings']['resultspath']
membranepdb = updated_data['simulation_parameters']['membrane_path_pdb']
topparpath = updated_data['simulation_parameters']['topparpath']  # toppar = topology + parameters
ligandsdict_path = updated_data['simulation_parameters']['ligandsdict_path'] 
modres_path = updated_data['simulation_parameters']['modres_path'] 

# Programs paths
chimera_path = detect_path("chimera")
vmd_path = detect_path("vmd")
psfgenpath = detect_path("psfgen")  # Download and use NAMDs psfgen for optimal results.
# Otherwise, you can always use the default HTMD one, but it has problems with organic halogens
os.environ['PATH'] += f":{chimera_path}:{vmd_path}:{psfgenpath}"

# Device machine settings
device_gpu = updated_data['machine_settings']['device_gpu']  # Which device GPU to use (0, 1, 2, 3...)
if not isinstance(device_gpu, int):
    raise Exception("The 'device_gpu' value must be an integer.")
acemd_path = "acemd"
acemd_queue = False
acemd_conda = f"source activate {scriptspath}htmd/"

# Preparation parameters
new_pdb_chain = updated_data['simulation_parameters']['preparation']['new_pdb_chain']
membrane_lipid_segid = updated_data['simulation_parameters']['preparation']['membrane_lipid_segid']
ph = updated_data['simulation_parameters']['preparation']['ph']
remove_sod2x50 = updated_data['simulation_parameters']['preparation']['remove_sod2x50']
coldist = updated_data['simulation_parameters']['preparation']['coldist'] # Distance bellow which two atoms are considered to collide
membrane_distance = updated_data['simulation_parameters']['preparation']['membrane_distance']# Distance between the pbc box and the space to be filled with membrane atoms
water_thickness = updated_data['simulation_parameters']['preparation']['water_thickness'] # Size in Z-axis of the solvation water layers
buffer = updated_data['simulation_parameters']['preparation']['buffer'] # Distance between solvation waters and the rest of the system
water_margin = updated_data['simulation_parameters']['preparation']['water_margin'] # Distance in the Z-axis to be penetrated by the solvation box 
                # to avoid the formation of a V O I D between the system and the solvation boxes
lipidlike_blacklist = updated_data['simulation_parameters']['preparation']['lipidlike_blacklist']    
detergent_blacklist = updated_data['simulation_parameters']['preparation']['detergent_blacklist']
glucids_blacklist = updated_data['simulation_parameters']['preparation']['glucids_blacklist']
blacklist = set(detergent_blacklist).union(set(glucids_blacklist))
tym_remove = updated_data['simulation_parameters']['preparation']['tym_remove']

# Topologies filenames and paths
toposfilenames = updated_data['simulation_parameters']['toposfilenames']

# Parameters filenames and paths
paramsfilenames = updated_data['simulation_parameters']['paramsfilesnames']

# Stream files (topology+parameters)
streamsfilenames = updated_data['simulation_parameters']['streamfilenames']

# Load topology, parameter, and stream files with our current basepath
topos = [os.path.join(topparpath, file) for file in toposfilenames]
params = [os.path.join(topparpath, file) for file in paramsfilenames]
streams = [os.path.join(topparpath, file) for file in streamsfilenames]

# Ligands for which parameters are not required (already exist in main CGenff)
AA={"LYR","ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
nucl = {"AMP","ADP","ATP","GMP","GDP","GTP","TMP","TDP","TTP","CMP","CDP","CTP"}
std_modres = {"SEP",'TPO'}
noparams_ligs = {'RET', 'CLR'}.union(AA).union(nucl).union(std_modres)

# Proteins containing this words in their name can be assumed to be GPCRs
gpcr_names = ['ceptor','rhodopsin','smoothened']
gpcr_uniprots = ['G1SGD4']

# Equillibration parameters
minim_steps = updated_data['simulation_parameters']['equilibration']['minim_steps']
equil_timestep = updated_data['simulation_parameters']['equilibration']['equil_timestep']
temperature = updated_data['simulation_parameters']['equilibration']['temperature']
equil_simtime = updated_data['simulation_parameters']['equilibration']['equil_simtime']
const_sel = updated_data['simulation_parameters']['equilibration']['const_sel']

# Simulation parameters
timestep = updated_data['simulation_parameters']['production']['timestep'] # Simulation timestep (femtoseconds)
trajperiod = updated_data['simulation_parameters']['production']['trajperiod'] # Simulation period: number of timesteps occuring between frames
repnum = updated_data['simulation_parameters']['production']['repnum'] # number of replicates 
prod_simtime = updated_data['simulation_parameters']['production']['prod_simtime'] # Production time per replicate

dummy_sel = 'name DUM'

# Load user input JSON file with system data
json_file=open(args.input_json)
json_str = json_file.read()
input_dict = json.loads(json_str)

# IMPORTANT: Peptide ligands must have L, L0, L1, L2, or PEP as its segid for
# this pipeline to work properly

print("Configuration and input files were successfully loaded!")