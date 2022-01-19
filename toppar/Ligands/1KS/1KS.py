
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/GPCR_simulations/Ligands/1KS/1KS.mol2") 
rc("addh") 
rc("setattr m name 1KS") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/GPCR_simulations/Ligands/1KS/1KS.mol2") 
            