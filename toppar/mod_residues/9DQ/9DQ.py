
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/GPCR_simulations/toppar/mod_residues/9DQ/9DQ.pdb") 
rc("delete element.H")
rc("addh") 
rc("setattr m name '9DQ'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/GPCR_simulations/toppar/mod_residues/9DQ/9DQ_chim.mol2") 
