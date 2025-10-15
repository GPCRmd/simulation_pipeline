import os
from config_pipeline import device_gpu, input_dict, resultspath, acemd_path, acemd_license, acemd_conda, acemd_queue
from simulate_structures_functions import define_production, job_commands, timestep, trajperiod, temperature, prod_simtime, repnum
from jobqueues.localqueue import LocalGPUQueue
#####################
## Part 4: Production
#####################

# Production protocol

# For each structure 
for entry in input_dict:    
    name = entry['name']
    pdbcode = entry['pdbcode']
    apo = entry['apo']
    # must match with equildir in equilibration launcher code and contain input and output of equilibration.
    modelname = name+'_apo' if apo else name
    equildir = resultspath+'equil/'+modelname+'/'

    if not os.path.exists(equildir):
        print("structure %s has not been yet equillibrated. Skipping...")
        continue
    for rep in range(1,repnum+1):
        
        try: 
            # If simulation for this PDB has already been run
            proddir = resultspath+'production/'+ modelname + '/rep_'+ rep +'/'
            if os.path.exists(proddir+'/output.xtc') or os.path.exists(proddir+'simrunning'):
                print("replicate %d of structure %s already has been simulated" %(rep, modelname))
                continue
            
            print('submitting replicate %d of %s' % (rep, modelname))
            if not os.path.exists(proddir):
                os.makedirs(proddir)
            
            # directory copy output of equilibration to production input (initial working directory for run_prod.sh).
            define_production(equildir, 
                              proddir, 
                              timestep, 
                              trajperiod, 
                              temperature, 
                              prod_simtime)

            # Prepare nohup command for ACEMD with conda environment activation
            command = (
            f"{acemd_path} --input {proddir} --device {device_gpu} > {proddir}nohup.out 2>&1"
            )
            with open(proddir + 'run.sh', 'w') as f:
                f.write('#!/bin/bash\n')
                f.write(f"{acemd_conda}\n")  # Activate conda environment
                f.write(command + '\n') # Run the acemd command
                f.write('exit')
            
            print(" ")    
            print(f"\033[94mStart running the equilibration script for {modelname}: {proddir}run.sh\033[0m")
            
            if not acemd_queue:
                print("Running equilibration...")
                # Execute the run.sh script
                os.system(f"{proddir}run.sh")
            else:
                print(f"ACEMD is set to local mode.")
                local = LocalGPUQueue()
                local.submit(proddir)
                local.wait()
            
        except Exception as e:
            print("model "+modelname+" could not be send to production because of ",e)
