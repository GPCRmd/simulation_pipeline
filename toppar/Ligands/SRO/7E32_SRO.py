
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/3th_simround//toppar/Ligands/SRO/7E32_SRO.pdb") 
rc("setattr m name 'SRO'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/3th_simround//toppar/Ligands/SRO/7E32_SRO.mol2") 
