import os
from config_pipeline import device_gpu, input_dict, resultspath, acemd_path, acemd_license, acemd_queue, acemd_conda
from simulate_structures_functions import define_equilibration, job_commands, const_sel, equil_simtime, minim_steps, equil_timestep, temperature
from jobqueues.localqueue import LocalGPUQueue

#########################
## Part 3: Equillibration
#########################

# Iterate by GPCRdb structures to simulate
for entry in input_dict:    
    name = entry['name']
    pdbcode = entry['pdbcode']
    apo = entry['apo']
    # Entry's data
    try:
        modelname = name+'_apo' if apo else name
        equildir = resultspath+'equil/'+modelname+'/'
        builddir = resultspath+'build/'+modelname+'/'
        if os.path.exists(equildir+'output.xtc') or os.path.exists(equildir+'simrunning'):
            print(" structure %s already has been equilibrated" % modelname)
            continue

        if not os.path.exists(equildir):
            os.makedirs(equildir)

        # Define equillibration parameters
        define_equilibration(
            builddir = builddir,
            outdir = equildir,
            parameters = builddir + 'parameters.prm',
            const_sel=const_sel,
            simtime = equil_simtime,
            minimize = minim_steps,
            timestep = equil_timestep, 
            temperature = temperature
        )
        
        # Prepare nohup command for ACEMD with conda environment activation
        nohup_command = (
        f"nohup {acemd_path} --input {equildir} --device {device_gpu} > {equildir}nohup.out 2>&1 &"
        )
        with open(equildir + 'run.sh', 'w') as f:
            f.write('#!/bin/bash\n')
            f.write(f"{acemd_conda}\n")  # Activate conda environment
            f.write(nohup_command + '\n') # Run the acemd command
        
        print(" ")    
        print(f"\033[94mStart running the equilibration script for {modelname}: {equildir}run.sh\033[0m")
        
        if not acemd_queue:
            # Make the script executable
            os.chmod(equildir + 'run.sh', 0o755)
            
            print("Running equilibration...")
            # Execute the run.sh script
            os.system(f"{equildir}run.sh")
        else:
            print(f"ACEMD is set to local mode.")
            local = LocalGPUQueue()
            local.submit(equildir)
            local.wait()
            
    except Exception as e:
        print("model "+modelname+" could not be send to equilibrate because of ",e)