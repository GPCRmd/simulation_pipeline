
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/GPCR_simulations/toppar/Ligands/7AB/7AB.pdb") 
rc("delete element.H")
rc("addh") 
rc("setattr m name 'COV'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/GPCR_simulations/toppar/Ligands/7AB/7AB_chim.mol2") 
