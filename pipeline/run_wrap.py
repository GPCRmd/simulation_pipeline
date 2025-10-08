import os
import multiprocessing as mp
from config_pipeline import input_dict, resultspath, strucpath, repnum
from simulate_structures_functions import transmem_atoms, wrap_alig_vmd

##########################
## Part 5: Wrap Structures
##########################

# Wrap trajectories obtained during production with an htmd command
prot_sel = "protein and segid P P0 P1 P2 P3 P4 P5 P6 P7 P8 P9" #Every chain in our system is assigned a PX segid. I dont think there will ever be more than 10
for entry in input_dict:    
    name = entry['name']
    pdbcode = entry['pdbcode']
    isgpcr = entry['isgpcr']
    prot_chain = entry['prot_chain']

    # Create a transmembrane selection (for alignment) 
    topopath = '%sproduction/%s/rep_1/structure.pdb' % (resultspath, modelname)
    transmem_sel = transmem_atoms(strucpath)
        
    modelname = name
    simpath = '%sproduction/%s/' % (resultspath, modelname)
    # Skip if all replicates are ready
    alldone=True
    for rep in range(1,repnum+1):
        trajpath = '%s/rep_%d/output_wrapped.xtc' % (simpath, rep)
        if not os.path.exists(trajpath):
            alldone=False
    if alldone:
        print("All replicates in dystem %s already wrapped. Skipping...."%(modelname))
        continue

    try:
        # Wrap individual replicates one by one
        pool = mp.Pool(3)
        for rep in range(1,repnum+1):
            
            trajpath = '%s/rep_%d/output.xtc' % (simpath, rep)
            if not os.path.exists(trajpath):
                print("System %s_%d not yet simulated. Skipping...."%(modelname,rep))
                continue
            
            outrajpath = '%s/rep_%d/output_wrapped' % (simpath, rep)
            if os.path.exists(outrajpath+'.xtc'):
                print("System %s_%d already wrapped. Skipping...."%(modelname,rep))
                continue
            
            # Wrap and align
            print("Wrapping and aligning %s..."% modelname)
            x = pool.apply_async(wrap_alig_vmd, args=(topopath, strucpath, trajpath, outrajpath, prot_sel, transmem_sel))
            print(x.get()) # Print errors

        pool.close()
        pool.join() 


    except Exception as e:
        print("System %s could not be wrapped because of %s" % (modelname,e))