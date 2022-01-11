
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/docking_project/toppar/Ligands/APE/ape_mod5_APE.pdb") 
rc("setattr m name 'APE'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/docking_project/toppar/Ligands/APE/ape_mod5_APE.mol2") 
