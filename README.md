# File description
+ **autosubmit**: automated script for the submission of multiple simulations into GPCRmd. Instructions contained inside
+ **demo**: Contains an example PDB file (4EJ4.pdb) to run the pipeline with, alongside a fulfilled inputs.json with its data. 
+ **expected_demo_output**: expected results after building and starting the equilibration of the demo set, included above.
+ **membrane/popc36_box_renumbered.pdb**: Patch of POPC membrane, necessary for the building of transmembrane systems.
+ **readme_afterbuild.txt**: Includes considerations to revise systems built with this pipeline, before its equilibration (e.g.: checking disulfide bridges, revising ligand molecules topology....).
+ **simulate_frompdb_structures.ipynb**: Main simulation pipeline for the building, equilibration, simulation and wrapping of systems. 
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
1. 
Downlaod files, put them in folder, and you are ready to go

# Instructions
+ Explain demo 

# Demo
+ Prepare input.json for 4EJ4, alonside structure
+ Describe what to do to run build pipeline
+ Add also compressed file with output of 4EJ4 build
+ Mention time to build usually (max. 15min)
+ Mention equillibration and proudcion times in our computer server
