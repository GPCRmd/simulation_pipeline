
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/GPCR_simulations/toppar/mod_residues/HRG/HRG.pdb") 
rc("delete element.H")
rc("addh") 
rc("setattr m name 'HRG'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/GPCR_simulations/toppar/mod_residues/HRG/HRG_chim.mol2") 
