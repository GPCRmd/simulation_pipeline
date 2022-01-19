
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/3th_simround//toppar/Ligands/5FW/7DHR_5FW.pdb") 
rc("setattr m name '5FW'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/3th_simround//toppar/Ligands/5FW/7DHR_5FW.mol2") 
