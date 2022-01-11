
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/3th_simround/toppar/Ligands/H98/6TKO_H98.pdb") 
rc("setattr m name 'H98'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/3th_simround/toppar/Ligands/H98/6TKO_H98.mol2") 
