
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/3th_simround/toppar/Ligands/PIO/6UP7_PIO.pdb") 
rc("setattr m name 'PIO'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/3th_simround/toppar/Ligands/PIO/6UP7_PIO.mol2") 
