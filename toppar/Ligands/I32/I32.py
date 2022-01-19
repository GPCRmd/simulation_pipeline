
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/GPCR_simulations/Ligands/I32/I32.mol2") 
rc("addh") 
rc("setattr m name I32") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/GPCR_simulations/Ligands/I32/I32.mol2") 
            