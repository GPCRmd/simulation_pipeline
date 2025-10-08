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
"""
import sys

# Passwords for paramchem
username = None
password = None
if not (username and password):
    raise Exception("Please define your password and username for CGenFF app on config_pipeline.py. It is required to parameterize small molecules")

# Paths to ACEMD3 and its license (REPLACE WITH YOUR OWN)
acemd_path = "acemd"
acemd_queue = True
acemd_conda = "source /home/agarcia/miniconda3/bin/activate htmd"
acemd_license = "SG_LICENSE_FILE=28000@tolkien.prib.upf.edu,ACELLERA_LICENSE_SERVER=28000@tolkien.prib.upf.edu"

import os

# Our main path
basepath = os.getcwd() + '/'
# Other Paths
scriptspath = basepath + 'pipeline/'
strucpath = basepath + 'input_structures/'
resultspath = basepath + 'simulation_output/'
membranepdb = basepath + 'membrane/popc36_box_renumbered.pdb'
topparpath = basepath + 'toppar/TOP_PARAMS_ACE3/'  # toppar = topology + parameters
ligandsdict_path = basepath + 'ligands.json'
modres_path = basepath + 'modified_residues.json'
slurmpath = basepath + 'fake_slurm/'

# Path to slurm queuing system binaries
# In our case, Ismael designed a bunch of small bash scripts (fake_slurm) which do ssh to Hydra and execute slurm there
# Remember to modify the "fake_slurm" commands by adding your username in them

# USER INPUTS: introduce a JSON dictionary with the data for your system in the following format
"""
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
"""

# Parameters
new_pdb_chain = 'R'
membrane_lipid_segid = 'MEM'
coldist = 1.3 # Distance bellow which two atoms are considered to collide
membrane_distance = 20 # Distance between the pbc box and the space to be filled with membrane atoms
water_thickness = 20 # Size in Z-axis of the solvation water layers
buffer = 2.4 # Distance between solvation waters and the rest of the system
water_margin = 4 # Distance in the Z-axis to be penetrated by the solvation box 
                 # to avoid the formation of a V O I D between the system and the solvation boxes
lipidlike_blacklist = {'OLA','OLB','OLC','PLM','HTG','LPP','PEF','2CV','SOG','TWT','STE','LMT',
                       'MPG','DGA','1WV','POV','FLC','PGW', 'PC1','LDA','J40','BNG'}
detergent_blacklist = {'OLA','OLB','PEG','GOL','BOG','OLC','P6G','P33','UNL','UNX','PLM','HTG',
                       '12P','LPP','PEF','2CV','SOG','TWT','PGE','SO4','STE','LMT','ACT','ACE',
                      'MHA','CIT','1PE','MPG','EPE','PG4','DGA','PO4','DMS','TAR','1WV','EDO',
                      'BU1','ACM','PG6','TLA','SCN','TCE','MES','EDT','POV','MLI','SIN','PGO',
                      'FLC','HTO','PGW','NO3', 'PE5', 'NH2', 'NME','NH4','IOD','ZN','BR','HG',
                      'PC1','LDA','C8E','J40','BNG'}
glucids_blacklist = {'MAN','NAG','BGC','TRE','9NS','BMA','FUC','A2G'}
blacklist = detergent_blacklist.union(glucids_blacklist)

# Ligands for which parameters are not required (already exist in main CGenff)
AA={"LYR","ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
nucl = {"AMP","ADP","ATP","GMP","GDP","GTP","TMP","TDP","TTP","CMP","CDP","CTP"}
std_modres = {"SEP",'TPO'}
noparams_ligs = {'RET', 'CLR'}.union(AA).union(nucl).union(std_modres)

# Proteins containing this words in their name can be assumed to be GPCRs
gpcr_names = ['ceptor','rhodopsin','smoothened']
gpcr_uniprots = ['G1SGD4']

# Equillibration parameters
minim_steps = 5000
equil_timestep = 2
temperature = 310
equil_simtime = 40
const_sel = 'protein and name C CA N O or not (protein or lipid or water or ions) and noh or segid IO WAT and noh'

# Simulation parameters
timestep = 4 # Simulation timestep (femtoseconds)
trajperiod = 50000 # Simulation period: number of timesteps occuring between frames
repnum = 3 # number of replicates 
prod_simtime = 500 # Production time per replicate

# Dummy atom stored in separated pdb
# It is used during removal of aromatic insertions by placing it on the middle of the ring and measuring distances  

dummy_sel = 'name DUM'

# Topologies filenames and paths
toposfilenames = [
                  'General_top_params/topologies/top_all36_cgenff.rtf',
                  'General_top_params/topologies/top_all36_prot.rtf',
                  'General_top_params/topologies/top_all36_na.rtf',
                  'General_top_params/topologies/top_all36_lipid.rtf',
                  'General_top_params/topologies/top_all36_carb.rtf',
                  'Specific_top_params/topologies/toppar_all36_lipid_ether.rtf',
                  'Specific_top_params/topologies/toppar_water_ions.rtf',
                  'Specific_top_params/topologies/toppar_all36_na_nad_ppi.rtf',
                  'Specific_top_params/topologies/toppar_all36_prot_na_combined.rtf',
                  'David_top_params/topologies/my_patches.rtf',
                  'David_top_params/topologies/toppar_all36_lipid_cholesterol_CLR.rtf',
                  'David_top_params/topologies/toppar_all36_prot_retinol.rtf',
                 ]

# Parameters filenames and paths
paramsfilenames = [
                   'General_top_params/parameters/par_all36_cgenff.prm',
                   'General_top_params/parameters/par_all36m_prot.prm',
                   'General_top_params/parameters/par_all36_na.prm',
                   'General_top_params/parameters/par_all36_lipid.prm',
                   'General_top_params/parameters/par_all36_carb.prm',
                   'Specific_top_params/parameters/toppar_all36_lipid_cholesterol.prm',
                   'Specific_top_params/parameters/toppar_all36_lipid_ether.prm',
                   'Specific_top_params/parameters/toppar_water_ions.prm',
                   'Specific_top_params/parameters/toppar_all36_na_nad_ppi.prm',
                   'Specific_top_params/parameters/toppar_all36_prot_na_combined.prm',
                   'Specific_top_params/parameters/toppar_all36_prot_retinol.prm',
                #    'David_top_params/parameters/fake_dihedrals.prm',
                   'David_top_params/parameters/modres_params.prm',
                   'David_top_params/parameters/modres_crossterm.prm',
                  ]

# Stream files (topology+parameters)
streamsfilenames = [
                ]
import json
json_file=open(basepath + "demo/inputs.json")
json_str = json_file.read()
input_dict = json.loads(json_str)
# IMPORTANT: Peptide ligands must have L, L0, L1, L2, or PEP as its segid for
# this pipeline to work properly


# Load topology, parameter, and stream files with our current basepath
topos = [os.path.join(topparpath, file) for file in toposfilenames]
params = [os.path.join(topparpath, file) for file in paramsfilenames]
streams = [os.path.join(topparpath, file) for file in streamsfilenames]

# Add Chimera's and VMD paths (replace by your own)
chimera_path = "/soft/system/software/Chimera/1.16/bin/"
vmd_path = "/soft/system/software/VMD/1.9.4a57/bin/vmd"
psfgenpath = None  # Download and use NAMDs psfgen for optimal results.
# Otherwise, you can always use the default HTMD one, but it has problems with organic halogens
os.environ['PATH'] += f":{chimera_path}:{vmd_path}:{slurmpath}"

