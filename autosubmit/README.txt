############################
## Manual for submit_sims.py
############################

This script automatically submits your simulation systems into GPCRmd. 
To do so, it requires an input json file with all the information necessary for a submission.
You can use "input_template.json" as template to create your own submission json.
This json must have an entry for each system you wish to submit, and each entry must have the following fields:

#### 1. Sys_info
Simulation specs and technical data, used in step1

	1. "name" : a name for your system.
	2. "type" : COMPLEX or APO depending on presence or absence of ligand in your simulated GPCR 
	3. "pdbid" : PDB ID of the structure from which your system is based. If your system was obtained by homology modelling, use HOMO here. If it comes from Alpha-fold, use ALPHA. Otherwise just use the closest-related existing PDB structure to your system.
	4. "description" : Short description of your simulation (optional. Leave blank if so wish).
	5. "source_type" : XRAY, NMR, EM (electron microscopy) or OTH (other). Source of simulated structure.  
	6. "id_dynamics_methods": simulation method used
		MM (molecular mechanics)
		QM (Quantum Mechanics)
		NMR 
		MONTECARLO (montecarlo simulation)
		MSM (MM MSM/HMM adaptative)
		METADYN (MM metadynamics)
		UMBRELLA (MM umbrella sampling)
		BIAS (MM biased)
		SCALED", 
	7. "software": Name of software used for simulation (ACMED3, GROMACS,...)
	8. "sversion": Version of the simulation software, 
	9. "ff": Force field used (e.g.: CHARMM)
	10. "ffversion": Version of the force field used 
	11. "id_assay_types": Type of simulation assay
		ORTO (orthosteric (un)binding)
		ACT (Activation/Inactivation)
		OLIGO (Oligomerization)
		ALLOBIND (Allosteric (un)binding)
		ALLOMOD (allosteric modulation)
		GPROTCOUP (G-protein (un)coupling)
		BARRCOUP (B-arrestin (un)coupling)
	12. "id_dynamics_membrane_types": Membrane type used
		IMP (Implicit)
		HOMO (Homogeneous)
		HETE (Heterogeneous)
	13. "id_dynamics_solvent_types": solvent type
		IMP (implicit)
		TIP3
		TIP4
		TIPS
		SPC
		SPCE (SPC/E)
		OTH (other)
	14. "timestep": simulation timestep
	15. "delta": simulation delta. Time lapse between frames in the simulation (in ns) 
	16. "add_info": Additional info (optional)
	
#### 2. Small molecules
Create an entry in this section for each non-protein molecule type in your system, including solvation waters lipids and ions. 
The format of these entries is the following
	1. "type" : role of this molecule in the simulation
		ORT (Orthosteric ligand)
		ALLO (allosteric ligand)
		EION (Experimental ions)
		ELIP (Experimental lipids)
		EWAT (Experimental waters)
		OTHC (other experimental compounds)
		BWAT (Bulk solvation waters)
		BLIP (Bulk lipids)
		BION (Bulk solvation lipids)
		BOTH" (Other bulk components) 
	2. "resname" : Residue name of this molecule in the simulated system
	3. "inchikey" : inchikey of this molecule.

### 3. Proteins
Create an entry for each protein in your system. By "protein" here we mean individual proteins, identifiable by Uniprots or equivlanets.
Ideally each individual protein should have its own chainID, but it is not required.
Should a protein chain be a chimera of two different proteins, make two entries and properly specify their coordinates.
	1. "isoform" : in case this protein is an isoform of the UNIPROT one, put here its number. Otherwise just put '1' 
	2. "uniprotkbac" : Uniprot ID of this protein Should your protein not have an uniprot leave this field blank
	3. "segs" : molecular coordinates of the residues belonging to this protein. Several entries can be created, in case protein is fragmented. 
		"chainid" : Chain ID as in simulation coordinates,
		"segid" : Segment IDs as in simulation coordinates,
		"beg" : Residue ID of the first residue in this segment,
		"end" : Residue ID of the last residue in this segment,
		"bond" : true/false, depending on wether or not this segment is bound to the previous one
		"source_type" : Experimental source of this segment structure. SHould be the same as in step1, except in the case of chimeras
			XRAY
			NMR 
			ABIN (ab initio modelling)
			HOMO (Homology modelling)
			THR (Threading)
			MD (Molecular Dynamics)
			EM (Electron microscopy)
			OTH (Other) 
	
	4. "isgpcr" : true/false wether this protein is or not a GPCR  
	5. "organism_id" : Uniprot ID of this protein's organism (unrequired if protein has uniprotkbac)
	6. "organism" : Uniprot name of this protein's organism (unrequired if protein has uniprotkbac)
	7. "name" : Name of the protein (unrequired if protein has uniprotkbac) 
	8. "other_names" : Other names this protein may have (separated by ';')(unrequired if protein has uniprotkbac)

### 4. Files
Files of your simulation.
	1. "coords" : full path to a  coordinates file of your simulated system, as in frame 1 of the simulation. 
	This file has must include all elements in your simulation (including solvation waters and membrane lipids).
	The proteins MUST have segment IDs identifying them, and preferibly chain IDs too. 	
	Accepted formats: PDB
	
	2. "topo" : full path to a topology file of your simulated system.
	Accepted formats: PSF, TOP

	3. "trajs" : list with full paths to the trajectory files of your simulation. Several files can be included here.
	Accepted formats: dcd, xtc, trj, trr
	
	4. "params" : Parameters used in your simulation. Stored for replicability
	Accepted formats: any, but must have '.prm' termination. If several files need to be uploaded, compress them in a 'tar.gz' file
	
	6. "other" : Any other files relevant to this simulaiton, including the protocol used to run it
	Accepted formats: any, as long as compressed in a 'zip' or 'tar.gz' file

### 5. References
Information of the publication in which this simulation will appear. 
This step can be omited, and manually fullfilled later on, if your submission is not yet published
	1. "doi" : a temporary placeholder can be included hear to bind toghether all simulations included in the same publication
	2. "authors" 
	3. "title" 
	4. "pmid" 
	5. "journal"
	6. "issue" 
	7. "volume" 
	8. "pages" 
	9. "publication_year" 
	10. "URL" 
