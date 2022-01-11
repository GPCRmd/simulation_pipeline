
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/PSYBIAS_simulations/toppar/Ligands/UN3/otava35_mut_UN3.pdb") 
rc("setattr m name 'UN3'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/PSYBIAS_simulations/toppar/Ligands/UN3/otava35_mut_UN3.mol2") 
