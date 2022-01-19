
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/3th_simround//toppar/Ligands/XGD/7L1V_XGD.pdb") 
rc("setattr m name 'XGD'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/3th_simround//toppar/Ligands/XGD/7L1V_XGD.mol2") 
