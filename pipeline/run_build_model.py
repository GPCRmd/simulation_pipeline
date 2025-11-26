from config_pipeline import input_dict, basepath, strucpath, resultspath, membranepdb, psfgenpath, topos, params, streams, strucpath, ph, tym_remove, remove_sod2x50
from simulate_structures_functions import gpcrdb_dict, blacklist, water_thickness, water_margin, buffer, coldist, membrane_distance, membrane_lipid_segid
from simulate_structures_functions import internal_waters, fix_and_prepare_input, make_apo, covalent_ligands
from simulate_structures_functions import get_opm, chimera_superimpose, prepare_system
from simulate_structures_functions import add_membrane, solvate_pdbmol, extra_parameters, get_caps
from simulate_structures_functions import renumber_resid_vmd, _recoverProtonations
from simulate_structures_functions import charmm
from htmd.ui import Molecule
import os
import numpy as np
import time
import traceback 
import json

###########################
## Part 2: Build the models
###########################

# Load modresdict and ligandsdict from JSON files
modresdict_path = os.path.join(strucpath, "modresdict.json")
ligandsdict_path = os.path.join(strucpath, "ligandsdict.json")

with open(modresdict_path, "r") as modres_file:
    modresdict = json.load(modres_file)

with open(ligandsdict_path, "r") as ligands_file:
    ligandsdict = json.load(ligands_file)

# Iterate by GPCRdb structures to simulate
for entry in input_dict:    
    try:

        # Entry's data
        name = entry['name']
        isgpcr = entry['isgpcr']
        pdbcode = entry['pdbcode']
        pdbfile = entry['pdbfile']
        curated = entry['curated']
        sod2x50 = entry['sod2x50']
        prot_chain = entry['prot_chain']
        gpcr_chain = entry['prot_chain'] if isgpcr else False
        apo = entry['apo']
        
        #Starting simulation
        start_time = time.time()        
        sysname = name+'_apo' if apo else name
        mystrucpath = strucpath+sysname+'/'
        os.makedirs(mystrucpath, exist_ok=True)

        # Skip if there is already a model build for this
        if os.path.exists(resultspath+'build/'+sysname+'/structure.pdb'):
            print('Build model for '+sysname+' already exists. Skipping...')
            continue

        # Check if simulation is aminergic
        if isgpcr:
            aminergic = gpcrdb_dict[pdbcode]['family'].startswith('001_001')
            adenosine = gpcrdb_dict[pdbcode]['family'].startswith('001_006_001')
        else:
            aminergic = None; adenosine = None
        
        # Add waters with homolwat if protein is a GPCR. Sodium 2x50 will also be added if a 
        # non-false pdbcode is added
        if isgpcr:
            (sod2x50, watered_filename) = internal_waters(pdbfile, pdbcode, gpcrdb_dict, apo, sod=sod2x50)
        else:
            watered_filename=pdbfile
        mol = Molecule(watered_filename)
        if remove_sod2x50: 
            # Remove any non-protein, non-ion, non-water thing on the system. 
            # Delete also sod2x50 (but we will keep 2x50 protonated)
            mol.remove('element Na')
        # Remove unnecessary ligand molecules: mostly crystalization detergents, quelants, buffers,
        # or post-traductional glicosilations
        mol.remove('resname '+' '.join(blacklist))
        
        # Remove 2x50Sodium from non-A-class GPCRs
        if isgpcr:
            if not gpcrdb_dict[pdbcode]['family'].startswith('001'):
                mol.remove('element NA')
                
        # Ismael's function to add labels (segid) for 'ligand' and 'protein' parts of the system
        mol_fixed,prot_segids = fix_and_prepare_input(mol,name,pdbcode,
                                                        modresdict,
                                                        isgpcr=isgpcr,
                                                        prot_chain=prot_chain,
                                                        tym_remove = tym_remove,
                                                        )
        # If the pipeline is running in 'apoform mode', remove any non-protein, non-ion, non-water thing on the system      
        # Delete also non-receptor proteins
        # If there's any, parameterize and rename covalent-bound ligands
        if apo:
            (mod_mol,prot_segids) = make_apo(mol_fixed,prot_chain)    
            covligs = []
        else:
            (mod_mol, covligs) = covalent_ligands(mol_fixed, name, ligandsdict)

        # Get aligned OPM structure
        thickness,opm_mol = get_opm(pdbcode)

        # Superimpose fixed molecule onto OPM counterpart for proper membrane fitting
        # Bit dirty, I know, but the best option avaliable
        modfile = 'mod_mol.pdb'
        mod_mol.write(modfile)
        opmfile = 'opm_mol.pdb'
        opm_mol.write(opmfile)
        aligfile = 'aligned_mol.pdb'
        chimera_superimpose(opmfile,modfile,aligfile)            
        mol_aligned = Molecule(aligfile)
        os.remove(modfile);os.remove(opmfile);os.remove(aligfile)
        
        #Center to receptor XY
        center = np.mean(mol_aligned.get('coords',sel='chain '+prot_chain),axis=0)
        mol_aligned.moveBy([-center[0],-center[1],0])
        
        # Prepare protein: asign titration states, flipping side chains of HIS, ASN and GLN; rotate some sidechains, optimize waters, etc.
        # Most of this is done with a HTMD function called proteinPrepare()
        # Skip step if we are working with curators structures
        prepared_mol = mol_aligned if curated else prepare_system(mol_aligned, pdbcode, thickness, gpcr_chain, sod2x50, aminergic, adenosine, ph)
        
        #Add membrane
        print('Adding membrane...')
        membranemol = Molecule(membranepdb)
        mol_membraned, membrane_resnames, membrane_segids, xreps, yreps = add_membrane(prepared_mol, membranemol,prot_segids,membrane_distance)

        #Solvate
        print('Solvating...')
        mol_solvated = solvate_pdbmol(mol_membraned,membrane_segids,water_thickness,water_margin,buffer=buffer,coldist=coldist,prefix='WT')

        #Obtain extra parameters for ligands and modified residues 
        ligstreams=extra_parameters(name, ligandsdict, modresdict, blacklist, covligs, basepath)
        
        # Assignign terminology for cap atoms of protein chain, depending if it is the receptor protein or not
        caps = get_caps(prot_segids, mol_solvated)
        #e.g.: {'P0': ['first ACE', 'last CT3'], 'P1': ['first ACE', 'last CT3']}

        #Pre-build model
        print('Pre-build...')
        prebuildmol = charmm.build(mol_solvated, 
                                    topo=topos, 
                                    param=params,
                                    stream=streams+ligstreams,
                                    caps=caps,
                                    outdir=resultspath+'/pre-build/'+sysname,
                                    ionize=False,
                                    psfgen=psfgenpath)

        # Save prebuild model topologies in files, and  store prebuild model in molecule object
        prebuild_psffile = prebuildmol.topoloc
        prebuild_pdbfile = os.path.splitext(prebuildmol.topoloc)[0]+'.pdb'
        prebuildmol = Molecule(prebuild_pdbfile, validateElements=False)
        _recoverProtonations(prebuildmol)

        # Checking of water/lipid ratio
        lipid_num = len(set(prebuildmol.get('resid',sel='segid '+membrane_lipid_segid)))
        solv_num = len(prebuildmol.get('index',sel='resname TIP3 and name OH2'))
        if float(solv_num) / lipid_num < 35:
            raise ValueError('Water/lipid ratio lower than 35.')

        #Renumber residues
        print('Renumbering...')
        mol_renumbered = renumber_resid_vmd(prebuildmol,'segid '+' '.join(membrane_segids),by=2)

        # Ionizing system
        print('Ionizing...')
        molbuilt = charmm.build(prebuildmol,
                                topo=topos, 
                                param=params,
                                stream=streams+ligstreams,                        
                                outdir=resultspath+'/ionize/'+sysname,
                                saltconc=0.15,
                                caps=caps,
                                psfgen=psfgenpath)

        build_psffile = molbuilt.topoloc
        build_pdbfile = os.path.splitext(molbuilt.topoloc)[0]+'.pdb'
        molbuilt = Molecule(build_pdbfile, validateElements=False)
        _recoverProtonations(molbuilt)

        #Building system
        print('Building...')
        molbuilt = renumber_resid_vmd(molbuilt,'segid "WT.*" or segid I',by=2)
        molbuilt = charmm.build(molbuilt, 
                                topo=topos, 
                                param=params,
                                stream=streams+ligstreams,                        
                                outdir=resultspath+'/build/'+sysname,
                                ionize=False,
                                psfgen=psfgenpath)

        print('End of %s after %s seconds\n' % (sysname, time.time() - start_time))

    except Exception as e:
        print("model "+sysname+" could not be build because ",e)
        print(traceback.format_exc())
