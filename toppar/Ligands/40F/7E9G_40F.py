
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/3th_simround/toppar/Ligands/40F/7E9G_40F.pdb") 
rc("setattr m name '40F'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/3th_simround/toppar/Ligands/40F/7E9G_40F.mol2") 
