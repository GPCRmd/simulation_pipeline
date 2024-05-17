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
DESCRIPTION INPUT JSON
AN EXAMPLE IS GIVEN IN DEMO

## Ligand parametrization
SOME MOLECULES PARAMETERS ARE NOT PRESENT IN CGENFF
WE USE PARAMCHEM TO OBTAIN ITS TOPLOGY AND PARAMETERS
PLEASE ENSURE LIGAND MOLECULES HAVE HYDROGENS FOR PROPER TOPOLOGY GUESSING. PARAMCHEM MIGHT FAIL OTEHRWISE

## BUILDING 
BUILDING OF SIMULABLE SYSTEMS FROM INPUT PDB
1.INTERNAL WATERS ADDED WITH HOMOLWAT (LINK)
2. REMOVE CRYSTALIZATION ARTIFACTS (BLACKLIST)
3. CHANGE RESIDUE AND SEGMENT NAMES
4. ALIGN STRUCTURE TO OPM COUNTERPART,FOR MEMBRANE ADDITION, USING CHIMERA
5. PREPARE SYSTEM WITH HTMD SYSTEMPREPARE: TRITATION STATES, SIDECHAIN ROTATION, PROTONATION STATES. SKIP IF PDB IS ALREADY PROPERLY PROTONATED
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
