import os
from config_pipeline import resultspath, input_dict
import MDAnalysis as mda
import MDAnalysis.analysis.rms as rms
import numpy as np

#####################
## EXTRA: ligand RMSD
#####################
### Assess quality of systems by calculating the RMSD of the ligand/s molecule/s

# Set paths and files
outpath = resultspath+'rmsd_lig_probe.tsv'
if os.path.exists(outpath):
    print('RMSDs already computed. SKipping...')
else:

    # Iterate over pdbcodes 
    results = []
    for entry in input_dict:    
        name = entry['name']
        pdbcode = entry['pdbcode']
        # Iterate over production replicates
        for mytrajid in ["1","2","3"]:
            
            try:
                sysname = name

                # Input files of simulation
                files_path = resultspath+'production/%s/rep_%s/'%(name, mytrajid)
                files_p = resultspath+'build/%s'%(name)
                mypdbpath = files_p+'structure.pdb'
                mypsfpath = files_path+'structure.psf'
                mytrajpath = files_path+'output_wrapped.xtc'

                # Skip if no trajectory
                if not os.path.exists(mytrajpath):
                    print("no trajectory replicate %s for system %s avalible. Skipping..."%(mytrajid, sysname))
                    continue

                print('computing ligand RMSD for trajectory %s of system %s' % (mytrajid, sysname))

                # Load trajectory and topology into MDA universe, and select protein atoms
                u = mda.Universe(mypsfpath, mytrajpath)
                ligsel = u.select_atoms("segid LIG and not resname CLR")

                # Compute rmsd, extract its values and put them in corresponding lists and dicts
                R = rms.RMSD(ligsel)
                R.run()
                rmsd = np.mean([ a[2] for a in R.rmsd ])
                results.append((sysname,rmsd)) 

            except Exception as e:
                print("error: system %s failed becasue %s"%(sysname,e))

    # Once everything is done, sort and write RMSD results
    out = open(outpath,'w')
    out.write("Simulated_system\tRMSD_ligand\n")
    results_sorted = sorted(results, key=lambda tup: tup[0])
    for line in results_sorted:
        out.write("%s\t%f\n"%(line[0],line[1])) 
    out.close()