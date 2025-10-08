import json
from config_pipeline import input_dict, basepath, username, password, scriptspath
from simulate_structures_functions import save_smalmol_mol2, get_lig_toppar, get_modres_toppar

################################
# Part 1: Ligand parametrization
################################

# Save mol2 files of ligand molecules and modres present in systems
(modresdict, ligandsdict, pdbfilesdict) = save_smalmol_mol2(input_dict, basepath, hydrogenate_ligands=False)

# Get topology-parameter files for ligands
get_lig_toppar(ligandsdict, basepath, username, password, pdbfiles=pdbfilesdict)

# Get topology-parameter files for modified residues
get_modres_toppar(modresdict, basepath, username, password, pdbfiles=pdbfilesdict)

# Save modresdict and ligandsdict to JSON files
with open(f"{scriptspath}/modresdict.json", "w") as modres_file:
    json.dump(modresdict, modres_file)

with open(f"{scriptspath}/ligandsdict.json", "w") as ligands_file:
    json.dump(ligandsdict, ligands_file)