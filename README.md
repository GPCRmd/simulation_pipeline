# File description
+ **autosubmit**: automated script for the submission of multiple simulations into GPCRmd. Instructions contained inside
+ **demo**: Contains an example PDB file (4EJ4.pdb) to run the pipeline with, alongside a fulfilled inputs.json with its data. 
+ **expected_demo_output**: expected results after building and starting the equilibration of the demo set, included above.
+ **membrane/popc36_box_renumbered.pdb**: Patch of POPC membrane, necessary for the building of transmembrane systems.
+ **readme_afterbuild.txt**: Includes considerations to revise systems built with this pipeline, before its equilibration (e.g.: checking disulfide bridges, revising ligand molecules topology....).
+ **simulate_frompdb_structures.ipynb**: Main simulation pipeline for the building, equilibration, simulation and wrapping of systems.
+ **simulate_structures_functions.py**: Functions, modules and variables used by the pipeline.
+ **toppar**: Contains a selection of CHARMM General Force Field parameters and topologies, necessary to run this pipeline (version c36 Jul 2021, publicly available at https://mackerell.umaryland.edu/charmm_ff.shtml).

# Software requirements
+ Linux-based OS (tested successfully in Ubuntu 20.04.6 LTS and CentOS Linux release 7.5.1804) 
+ Visual Molecular Dynamics (VMD, must be available from the PATH)
+ UCSF Chimera (again, must be available from PATH)
+ Python (min. 3.10) and the following modules:
  + HTMD
  + MDAnalyisis
  + Biopython
  + Beautiful Soup
  + requests

# Installation and usage instructions
1. Download the content of this repo, and extract its contents into a folder.
2. Install the software and python modules specified above.
3. Open "simulate_frompdb_structures.ipynb as a jupyter notebook

## Input preparation
This pipeline uses a .json file as main input. This json must have the following structure:

```
[
    {
        "name" : "whatever name you wish for your system",
        "pdbfile" : "/path/to/your/input/structure.pdb",
        "modres" : ["residue name 1","residue name 1"], (all residue names for any non-cannonical protein residues present in your sequence) 
        "ligands" : [{
                        "resname":"residue name of your non-protein ligand molecule in the PDB file",
                        "name" : "name of this molecule",
                        "covalently_bound" : false/True (is the ligand covalently bound to the protein???), 
                        "inchikey" : "inchikey of this molecule"
                    },
                    {
                        "resname":"residue name of your second non-protein ligand molecule in the PDB file",
                        "name" : "name of this molecule",
                        "covalently_bound" : false/True (is the ligand covalently bound to the protein???), 
                        "inchikey" : "inchikey of this molecule"
                    },
                    .....],
        "apo": false/true (do you wish to simulate an apo-version of this structure, by removing all ligands and not-main-proteins of it??),
        "prot_chain" : "chain id of the main protein in this system (a GPCR, for us)",
        "pdbcode" : "pdb ID of the structure you are trying to run (or its most similar counterpart)", 
        "curated": false/true (if false, the structure will be reprotonated during the building process. Set as true if you wish to preserve your manual protonations), 
        "sod2x50": true/false (do you wish to add the conserved sodium near 2x50??(added using homolwat)), 
        "isgpcr": true/false (does this system contain a GPCR?) 
    },
    {
      (Add as many systems as you wish to build and simulate)
    },........
]
```
An input.json example is included in the demo folder

## Force field and ligand parametrization
To build and simulate our systems we always use the most recent version of CHARMM General Force Field (CGenFF, https://mackerell.umaryland.edu/charmm_ff.shtml#charmm). A partial copy of this force field is included in 'toppar/' to run the demo set. 
It is most likely that your ligand molecules are not included in CGenFF. For that reason, this pipeline automatically obtains both parameters and topologies for all specified ligand molecules in your structure using "paramchem" (https://cgenff.silcsbio.com/initguess/). The resulting '.str' file for each ligand will be saved in toppar/Ligands/
**WARNING**: ligand molecules in your PDB structure **must** be properly protonated and hydrogenated. Otherwise, parametrization is likely to fail, returning an empty '.str' file.

## System building
The third cell in 'simulate_frompdb_structures.ipynb' builds a system for every entry in input.json. The main tasks performed in this part of the script are:
1. Addition of internal GPCR waters (if system is a GPCR) through homology using HomolWat web server ([LINK](https://alf06.uab.es/homolwat/))
2. Remove potentially unwanted co-crystallized molecules. This mainly adresses detergents and other compoonents required for crystallization, but of no biological significance. The exact list of 'blacklisted' molecules we remove is included in 'simulate_structures_functions.py'.
3. Setting of segment IDs for each component in the system, and some residue and atom renaming for the proper funcitoning of charmm.build() in successive steps
4. Alignment of your input structure its equivalent in Orientations of Proteins in Membranes database (OPM, https://opm.phar.umich.edu/) using USCF chimera. This is required for proper membrane placing and requires the PDB id of an at least somewhat similar structure present in OPM.
5. If "curated" is 'false' we apply HTMD's systemprepare function (https://software.acellera.com/htmd/tutorials/protein-preparation.html). This function adresses matters such as assigning titration states at the user-chosen pH; flipping the side chains of HIS, ASN, and GLN residues; and optimizing the overall hydrogen bonding network. 
6. ADD MEMBRANE FROM EXISTING MEMBRANE PATCH
7. SOLVATE IN A WATER BOX. PARAMETERS OF WATER BOX INCLUDED IN SIMULATE_STRUCTURES_FUNCTIONS.PY
8. ADD CAPS: ONLY TO MAIN PROTEIN (PEPTIDE LIGANDS EXCLUDED)
9. 3 BUILD ROUNDS WITH CHARMM.BUILD() : FIRST TO OBTAIN THE SYSTEM TOPOLOGY, SECOND TO IONIZE SYSTEM, THIRD TO ACTUALLY BUILD THE SYSTEM. 
EACH BUILD SAVED SEPARATEDLY IN SIMULATION_OUTPUT

## equiliBRATION
DEFAULT PARAMETERS IN SIMULATE_STRUCTURES_FUNCTIONS()
WE USE ACEMD3 BY DEFAULT AS SIMULATION SOFTWARE
SYSTEMS WILL BE PREPARED TO RUN A EQUILLIBRATION IN ACEMD3. USE RUN.SH SCRIPT
EQUILLIBRATIONS AND SIMULATIONS ARE EXPECTED TO RUN IN CLUSTER, NOT WORKSTATIONS. 
PROVIDED DEMO TOOK 1-2 DAYS OF COMPUTATION TO BE EQUILLIBRATIED IN CLUSTER, WITH A GPU GTX 1080Ti
ONce run is done, files will be in simulation_output

## Production
ACTUAL SIMULATION STEP
3 REPLICATES BY DEFAULT, CHANGE AT WILL
AGAIN, EXPECTED TO RUN IN CLUSTER WITH ACEMD3. USE RUN.SH SCRIPT
IN OUR CLUSTER IT TOOK 5 DAYS

## Wrapping and aliging of simulated systems
Depends on VMD
Aligns systems around transmembrane protein atoms

## Other
A ligand RMSD calculator is provided, to ensure stability of ligand
Another script is provided to submit your simulation into GPCRmd, if necessary. A password and username must be specified first


# Instructions
+ Explain demo 

# Demo
+ Prepare input.json for 4EJ4, alonside structure
+ Describe what to do to run build pipeline
+ Add also compressed file with output of 4EJ4 build
+ Mention time to build usually (max. 15min)
+ Mention equillibration and proudcion times in our computer server
