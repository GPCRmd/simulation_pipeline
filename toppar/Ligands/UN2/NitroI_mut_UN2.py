
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/PSYBIAS_simulations/toppar/Ligands/UN2/NitroI_mut_UN2.pdb") 
rc("setattr m name 'UN2'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/PSYBIAS_simulations/toppar/Ligands/UN2/NitroI_mut_UN2.mol2") 
