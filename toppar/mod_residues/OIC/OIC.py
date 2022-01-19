
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/GPCR_simulations/toppar/mod_residues/OIC/OIC.pdb") 
rc("delete element.H")
rc("addh") 
rc("setattr m name 'OIC'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/GPCR_simulations/toppar/mod_residues/OIC/OIC_chim.mol2") 
