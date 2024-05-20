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
IMPORTANT: if your system includes a peptide ligand, please set "L" as its chain ID in the input PDB file.

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
6. Membrane addition: since this system has been already aligned to membrane orientation, the membrane is now creating by adding tiles of our defined membrane patch (in 'membrane/popc36_box_renumbered.pdb') and removing those lipids coliding with the protein. The membrane patch we provide is POPC-only, feel free to change it by any membrane composition of your preference. Parameters (included in 'simulate_structures_functions.py'):
  + _coldist_: minimum distance at which two atoms are considered to colide
  + _membrane_distance_: 'size' of the membrane, understood as distance from the protein to the edge of the box (in x/y direction)
8. Solvate system, including it in a TIP3 waterbox. Parameters:
  + _water_thickness_: Distance in the Z-axis from the membrane to the top/bottom of the PBC box
  + _water_margin_: Distance in the Z-axis to be penetrated by the solvation box
  + _buffer_: Distance between solvation waters and the rest of the system
9. Addition of cappings to the system's protein chains. By default we use ACE and CT3 cappings for the main proteins, and NTER and CTER (=no capping) for peptide ligands. **IMPORTANT**: Any protein sequence with chain ID L is assumed to be a peptide ligand.
10. Actual "building" of the system using charmm.build() function integrated in HTMD and CGenFF topologies. Systems are build in three rounds:
    1. A first one to stablish the system topology
    2. A second, to properly ionize the build system
    3. A third, to fully build the system

Built systems are saved in simulation_output/build/name_of_your_system

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
