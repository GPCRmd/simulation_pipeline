
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/3th_simround//toppar/Ligands/FFR/7C2E_FFR.pdb") 
rc("setattr m name 'FFR'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/3th_simround//toppar/Ligands/FFR/7C2E_FFR.mol2") 
