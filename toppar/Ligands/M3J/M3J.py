
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open /gpcr/users/daranda/doctorat/GPCR_simulations/Ligands/M3J/M3J.mol2") 
rc("addh") 
rc("setattr m name M3J") 
writeMol2(openModels.list(modelTypes=[Molecule]), "/gpcr/users/daranda/doctorat/GPCR_simulations/Ligands/M3J/M3J.mol2") 
            