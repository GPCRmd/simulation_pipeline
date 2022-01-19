
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/3th_simround/toppar/mod_residues/TYS/TYS.pdb") 
rc("delete element.H")
rc("addh") 
rc("setattr m name 'TYS'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/3th_simround/toppar/mod_residues/TYS/TYS_chim.mol2") 
