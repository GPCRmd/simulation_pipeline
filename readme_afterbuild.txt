#####################
### Curation protocol
#####################

This pipeline is a very nice tool, yet it presents limitations. In this document, we present a list of features that, ideally, should be revised and checked manually in systems build with this pipeline before simulating them:

###  1. INTERNAL WATER PLACEMENT
Performed by using the online web server HomolWat in the pipeline. 
 - Are there any clashing waters?
 - Are there any missing waters required for this receptor in particular?

###  2. SODIUM ION AT ALLOSTERIC SITE D2.50
The 2.50 allosteric sodium is placed using HomolWat too. However, this software is not enterily reliable and sometimes is unable to place it, specially for newer structures.
 - In case you requested it, is there a properly placed sodium ion near 2.50?

### 3. PROTONATION/TAUTOMERIC STATES
Protonation and tautomeric states have been automatically assigned by PROPKA in the pipeline (througth HTMD's "proteinPrepare()" function).
 - The pipeline should have protonated esidue D2.50 in absence of sodium ion, and de-protonate it in its presence. Is it correct?
 - PROPKA might not protonate properly residues facing the lipid bilayer. Ensure its proper protonation manually if you so wish.

### 4. DISULFIDE BONDS
Througth the "proteinPrepare()" function, our pipeline assigns the resname CYX to any two cysteines close enougth to form a disulfide bridge.
In the charmm.build() phase, disulfide bridges will be assigned to any pair of CYS close enougth or to pairs of CYX regardless of distance.
 - Are all the disulfide bridges required in your receptor properly set? If not, you might need to replace them by CYX in the initial structure

### 5. LIGAND PROTONATION AND PARAMETRIZATION
"Ligand" and "modres" molecules are usually not fitted for PROPKA, so the pipeline uses chimera to protonate them and paramchem to obtain its parameters and topology. 
However, this approach is rather imperfect, and chimera often gets the topology of non-standard molecules wrong.
 - We encourage you to manually place hydrogens for this kind of molecules, and if you choose to do it automatically please check the result.
   