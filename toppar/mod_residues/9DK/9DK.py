
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/GPCR_simulations/toppar/mod_residues/9DK/9DK.pdb") 
rc("delete element.H")
rc("addh") 
rc("setattr m name '9DK'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/GPCR_simulations/toppar/mod_residues/9DK/9DK_chim.mol2") 
