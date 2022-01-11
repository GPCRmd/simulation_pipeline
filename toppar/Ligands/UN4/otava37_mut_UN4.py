
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/PSYBIAS_simulations/toppar/Ligands/UN4/otava37_mut_UN4.pdb") 
rc("setattr m name 'UN4'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/PSYBIAS_simulations/toppar/Ligands/UN4/otava37_mut_UN4.mol2") 
