import json
import os
import traceback
import re
import requests
import zipfile,io
import glob
import numpy as np
from bs4 import BeautifulSoup   
import random
import tempfile
import time
# import MDAnalysis as mda
# from MDAnalysis.analysis import rms
from Bio import pairwise2
from Bio.Align import substitution_matrices 
import tarfile
import multiprocessing as mp
from requests_toolbelt.multipart.encoder import MultipartEncoder

#from requests_toolbelt import MultipartEncoder
from subprocess import Popen, PIPE

#HTMD things
os.environ["HTMD_NONINTERACTIVE"] = "1"
import htmd
from htmd.ui import *
from moleculekit.tools.sequencestructuralalignment import sequenceStructureAlignment
try:
    from htmd.protocols.equilibration_v3 import Equilibration
except Exception:
    from htmd.protocols.equilibration_v2 import Equilibration
from htmd.protocols.production_v6 import Production
from htmd.builder.builder import removeLipidsInProtein, tileMembrane, minimalRotation,removeAtomsInHull
from moleculekit.util import rotationMatrix, sequenceID
# from moleculekit.opm import get_opm_pdb
from htmd.config import config
from htmd.builder.charmm import _recoverProtonations

# IMPORTANT!!!! vmd needs to be installed for this pipeline to run properly
# config(viewer='vmd')

####################
## Initial variables
####################

# Get current GPCRdb data into a Json
gpcrdb_data = requests.get('http://gpcrdb.org/services/structure/').json()
gpcrdb_dict = { entry['pdb_code'] : entry for entry in gpcrdb_data }

# Lists of GPCR pdbs
first_round = {'1u19', '2rh1', '2y00', '2y02', '2y03', '2y04', '2ycw', '3d4s', '3dqb', '3ny8', '3ny9', '3nya', '3odu', '3oe0', '3pbl', '3pds', '3pqr', '3rze', '3uon', '3v2y', '3vw7', '3zpq', '3zpr', '4ami', '4amj', '4bvn', '4djh', '4dkl', '4ea3', '4grv', '4iaq', '4iar', '4ib4', '4k5y', '4l6r', '4lde', '4ldl', '4ldo', '4mbs', '4mqs', '4mqt', '4n6h', '4oo9', '4or2', '4phu', '4pxz', '4py0', '4qkx', '4rwd', '4rws', '4s0v', '4u15', '4u16', '4xee', '4xnv', '4xnw', '4xt1', '4yay', '4z34', '4z35', '4z36', '4zj8', '4zjc', '4zud', '5a8e', '5c1m', '5cgc', '5cgd', '5cxv', '5dhg', '5dhh', '5dsg', '5glh', '5jqh', '5l7d', '5l7i', '5tgz', '5u09', '5uen'}
second_round = {"3C9L","3C9M","6ZDR","6ZDV",'6RZ8', '6FFI', '6KJV', '6TQ6', '6DRY', '5WIV', '5V57', '6LI2', '6MH8', '5MZJ', '5KW2', '5X33', '5ZHP', '6PS6', '6FJ3', '6D35', '4NTJ', '6A94', '6AK3', '5UVI', '6FFH', '5WIU', '6KQI', '5XR8', '6RZ9', '5VEX', '6LW5', '6FKA', '6PS7', '6M9T', '6KNM', '6ME6', '6GPX', '5XSZ', '5K2A', '5YQZ', '6PS0', '5K2C', '6ME7', '2Z73', '6D32', '6RZ4', '5YC8', '5T1A', '6AKY', '4N4W', '5WS3', '6LRY', '6J21', '5TE5', '4O9R', '6FK6', '1GZM', '6DRX', '6TQ7', '5NDZ', '6IGL', '5VRA', '5ZBQ', '6TPN', '6DS0', '5NLX', '5X93', '6HLP', '5OLG', '5OM4', '5O9H', '5XPR', '6A93', '5GLI', '6OBA', '6D27', '6HLO', '6RZ5', '5TE3', '6IIV', '6KP6', '6TQ4', '6TPJ', '5V56', '5UNG', '6FK9', '6PS8', '2YCZ', '6JZH', '5XRA', '5ZKB', '4JKV', '5MZP', '5K2D', '5OLV', '5X7D', '5T04', '5OLO', '5DYS', '6FK7', '6J20', '6K1Q', '5EE7', '6D26', '6PS1', '6DRZ', '4GBR', '6FKC', '6TP6', '6LUQ', '5ZKQ', '5UNH', '6PRZ', '5UNF', '5NX2', '6LI0', '5NM2', '5D6L', '6AQF', '6KK1', '6C1Q', '4Z9G', '6LI1', '5WF5', '6FKB',  '5OLH', '4XT3', '5ZK3', '5VEW', '4XES', '6HLL', '6TOT', '5N2R', '5WF6', '5K2B', '6TPG', '6PS4', '5F8U', '6ME4', '6GT3', '6RZ6', '6ME8', '6ME2', '5VBL', '5N2S', '5ZBH', '6PS3', '6IIU', '4NC3', '5JTB', '5ZKP', '6TQ9', '6GPS', '6TP3', '6RZ7', '6ME3', '5ZK8', '6C1R', '5V54', '5UIG', '5TVN', '5OM1', '6AKX', '5ZKC', '5NM4', '5TZY', '6TOD', '4GPO', '6TO7', '6CM4', '6TOS', '6PS5', '6ME9', '5NDD', '6KPC', '6KK7', '4BUO', '6OL9', '6PS2', '6FK8', '6ME5', '6TP4', '4QIM', '6IGK', '5WQC', '5UIW', '5TZR', '5OLZ', '4EJ4', '6E59'}

# Parameters
new_pdb_chain = 'R'
membrane_lipid_segid = 'MEM'
coldist = 1.3 # Distance bellow which two atoms are considered to collide
membrane_distance = 20 # Distance between the pbc box and the space to be filled with membrane atoms
water_thickness = 20 # Size in Z-axis of the solvation water layers
buffer = 2.4 # Distance between solvation waters and the rest of the system
water_margin = 4 # Distance in the Z-axis to be penetrated by the solvation box 
                 # to avoid the formation of a V O I D between the system and the solvation boxes
lipidlike_blacklist = {'OLA','OLB','OLC','PLM','HTG','LPP','PEF','2CV','SOG','TWT','STE','LMT',
                       'MPG','DGA','1WV','POV','FLC','PGW', 'PC1','LDA','J40','BNG'}
detergent_blacklist = {'OLA','OLB','PEG','GOL','BOG','OLC','P6G','P33','UNL','UNX','PLM','HTG',
                       '12P','LPP','PEF','2CV','SOG','TWT','PGE','SO4','STE','LMT','ACT','ACE',
                      'MHA','CIT','1PE','MPG','EPE','PG4','DGA','PO4','DMS','TAR','1WV','EDO',
                      'BU1','ACM','PG6','TLA','SCN','TCE','MES','EDT','POV','MLI','SIN','PGO',
                      'FLC','HTO','PGW','NO3', 'PE5', 'NH2', 'NME','NH4','IOD','ZN','BR','HG',
                      'PC1','LDA','C8E','J40','BNG'}
glucids_blacklist = {'MAN','NAG','BGC','TRE','9NS','BMA','FUC','A2G'}
blacklist = detergent_blacklist.union(glucids_blacklist)

# Ligands for which parameters are not required (already exist in main CGenff)
AA={"LYR","ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
nucl = {"AMP","ADP","ATP","GMP","GDP","GTP","TMP","TDP","TTP","CMP","CDP","CTP"}
std_modres = {"SEP",'TPO'}
noparams_ligs = {'RET', 'CLR'}.union(AA).union(nucl).union(std_modres)

# Proteins containing this words in their name can be assumed to be GPCRs
gpcr_names = ['ceptor','rhodopsin','smoothened']
gpcr_uniprots = ['G1SGD4']

# Equillibration parameters
minim_steps = 5000
equil_timestep = 2
temperature = 310
equil_simtime = 40
const_sel = 'protein and name C CA N O or not (protein or lipid or water or ions) and noh or segid IO WAT and noh'

# Simulation parameters
timestep = 4 # Simulation timestep (femtoseconds)
trajperiod = 50000 # Simulation period: number of timesteps occuring between frames
repnum = 3 # number of replicates 
prod_simtime = 500 # Production time per replicate

# Dummy atom stored in separated pdb
# It is used during removal of aromatic insertions by placing it on the middle of the ring and measuring distances  

dummy_sel = 'name DUM'

# Topologies filenames and paths
toposfilenames = [
                  'General_top_params/topologies/top_all36_cgenff.rtf',
                  'General_top_params/topologies/top_all36_prot.rtf',
                  'General_top_params/topologies/top_all36_na.rtf',
                  'General_top_params/topologies/top_all36_lipid.rtf',
                  'General_top_params/topologies/top_all36_carb.rtf',
                  'Specific_top_params/topologies/toppar_all36_lipid_ether.rtf',
                  'Specific_top_params/topologies/toppar_water_ions.rtf',
                  'Specific_top_params/topologies/toppar_all36_na_nad_ppi.rtf',
                  'Specific_top_params/topologies/toppar_all36_prot_na_combined.rtf',
                  'David_top_params/topologies/my_patches.rtf',
                  'David_top_params/topologies/toppar_all36_lipid_cholesterol_CLR.rtf',
                  'David_top_params/topologies/toppar_all36_prot_retinol.rtf',
                 ]

# Parameters filenames and paths
paramsfilenames = [
                   'General_top_params/parameters/par_all36_cgenff.prm',
                   'General_top_params/parameters/par_all36m_prot.prm',
                   'General_top_params/parameters/par_all36_na.prm',
                   'General_top_params/parameters/par_all36_lipid.prm',
                   'General_top_params/parameters/par_all36_carb.prm',
                   'Specific_top_params/parameters/toppar_all36_lipid_cholesterol.prm',
                   'Specific_top_params/parameters/toppar_all36_lipid_ether.prm',
                   'Specific_top_params/parameters/toppar_water_ions.prm',
                   'Specific_top_params/parameters/toppar_all36_na_nad_ppi.prm',
                   'Specific_top_params/parameters/toppar_all36_prot_na_combined.prm',
                   'Specific_top_params/parameters/toppar_all36_prot_retinol.prm',
                #    'David_top_params/parameters/fake_dihedrals.prm',
                   'David_top_params/parameters/modres_params.prm',
                   'David_top_params/parameters/modres_crossterm.prm',
                  ]

# Stream files (topology+parameters)
streamsfilenames = [
                ]

# Url for submissions
mainurl = 'https://www.gpcrmd.org'

# Threshold activation distances (obtained from analyzing all GPCRdb refined structures)
#with open('./threshold_dists.pickle', 'rb') as handle:
#    activation_threholds = pickle.load(handle)


####################
## David's functions
####################

def load_mol_pdbcode(pdbcode):
    """
    Download a PDB id and load it into a Molecule object
    """

    url = 'https://files.rcsb.org/download/%s.pdb' % pdbcode
    pdb_data = requests.get(url).content

    # Save to temporary file if needed, or load via Molecule
    with open('temp.pdb', 'wb') as f:
        f.write(pdb_data)

    mol = Molecule('temp.pdb')
    os.system('rm temp.pdb')

    return(mol)

def ligands_by_system(ligandsdict):
    """
    Creates a json with all the systems codes and their ligand molecules 
    """
    with open('ligands_by_system_new.json', 'w') as ou:
        dasdict = {}
        for pdb in ligandsdict:
            dasdict[pdb] = list(ligandsdict[pdb].keys())
        json.dump(dasdict, ou, indent= 4)

def json_dict(path):
    """
    Converts json file to pyhton dict.
    """
    json_file=open(path)
    json_str = json_file.read()
    json_data = json.loads(json_str)
    return json_data

def get_complexsignal_pdbs():
    """
    Obtain PDB codes of GPCRdb refined structures with coupled auxiliary proteins
    Also get the name of this GPCR family
    """
    
    # Download crude HTML data from gpcrdb/structure and process its table into a list of lists
    html_data = requests.get('http://gpcrdb.org/structure/').text
    table_data = [[cell.text.replace('\n','') for cell in row("td")]
                             for row in BeautifulSoup(html_data)("tr")]

    # From this list of lists, extract information about the signaling proteins associated to each structure
    gprdb_extradata = { entry[7] : { 'signal_prot' : entry[14], 'family' : entry[3], 'type' : entry[2] } for entry in table_data 
                       if (len(entry) >= 8) and (len(entry[7]) == 4) }

    # Get PDB codes of structures with signaling proteins
    complex_signal_pdbset = { pdbcode for pdbcode,data in gprdb_extradata.items() if data['signal_prot'] != '-' }

    # Classify obtained GPCRs pdbcodes by family and type
    complexsignal_byfamily = {}
    for pdbcode,data in gprdb_extradata.items():
        family = data['family']
        gpcr_type = data['type'] 
        signal = data['signal_prot']

        # Skip those not bound to a signal protein
        if signal == '-':
            continue
        try:
            complexsignal_byfamily[family][gpcr_type].append(pdbcode)
        except KeyError:
            try:
                complexsignal_byfamily[family][gpcr_type] = [pdbcode]
            except KeyError:
                complexsignal_byfamily[family] = { gpcr_type : [pdbcode] }
            

    return(gprdb_extradata, complex_signal_pdbset, complexsignal_byfamily)

def get_GPCRdb_nonsimulated(gpcrdb_dict):
    """
    Returns a list of PDB codes from the GPCRdb refined structures not yet simulated in GPCRmd
    """
    
    #Make set with the pdb codes of the structures in GPCRdb
    gpcrdb_pdbs = set(gpcrdb_dict.keys())

    # Load a random name-to-dyn json from contactmaps
    # This Jsons contain a dictionary with the dynIDs and the full names of the GPCR simulated
    # This way I can get all the pdb codes of the GPCRs presents in GPCRmd
    response = requests.get('http://www.gpcrmd.org/dynadb/files/Precomputed/get_contacts_files/contmaps_inputs/all/cmpl/lg/name_to_dyn_dict.json')
    soup = BeautifulSoup(response.text, 'html.parser')
    gpcrmd_sims = json.loads(str(soup))

    # Take the pdb code from each full name in the json, and store in set
    gpcrmd_pdbs = set()
    pdb_pat = re.compile('\((\w+).*\) \(.*\)$')# Take objects whatever is inside of the first parenthesis
    for sim in gpcrmd_sims:
        match_pdb = re.search(pdb_pat, sim[1])
        if match_pdb:
            gpcrmd_pdbs.add(match_pdb.group(1))

    # Get to-be-simulated GPCR pdbs. That is the ones that are in GPCRdb but not in GPCRmd
    not_simulated = gpcrdb_pdbs - gpcrmd_pdbs

    return not_simulated

def auld_download_GPCRdb_structures(pdb_set, strucpath):
    """
    Download (if they exist) the refined GPCRdb structures for the pdb codes in the pdb_set.
    PDB codes without a refined structure will be removed from pdb_set
    """
    pdb_set_nonrefined = set()
    set_length = len(pdb_set)
    i = 0
    for pdbcode in pdb_set:
        mystrucpath = strucpath+pdbcode+'/'
        os.makedirs(mystrucpath, exist_ok = True)
        i += 1
        
        print('Downloading %s structure (%d/%d)' % (pdbcode, i, set_length))
        # If files for this simulation already exists
        if glob(mystrucpath+'*pdb'):
            print('Structure for %s already present. Skipping...' % pdbcode)
        else:
            # Try two possible GPCRdb repositories to download the refined structure
            response = requests.get('https://gpcrdb.org/structure/homology_models/'+pdbcode+'_refined_full/download_pdb', stream=True)
            if not response.ok:
                response = requests.get('http://build.gpcrdb.org/structure/homology_models/'+pdbcode+'_refined_full/download_pdb', stream=True)
            if not response.ok:
                pdb_set_nonrefined.add(pdbcode)
                print('could not download %s refined structure. Skipping...' % (pdbcode))
            else:
                #Extract downloaded files
                zippy = zipfile.ZipFile(io.BytesIO(response.content))
                zippy.extractall(mystrucpath)

def download_GPCRdb_structures(pdb_set, strucpath):
    """
    Download (if they exist) the refined GPCRdb structures for the pdb codes in the pdb_set.
    PDB codes without a refined structure will be removed from pdb_set
    """
    pdb_set_nonrefined = set()
    set_length = len(pdb_set)
    refined_baselink = 'https://gpcrdb.org/structure/refined/'
    i = 0
    for pdbcode in pdb_set:
        mystrucpath = strucpath+pdbcode+'/'
        os.makedirs(mystrucpath, exist_ok = True)
        i += 1
        
        try:
            print('Downloading %s structure (%d/%d)' % (pdbcode, i, set_length))
            # If files for this simulation already exists
            outpdb = mystrucpath+pdbcode+'_refined.pdb'
            if os.path.exists(outpdb):
                print('Structure for %s already present. Skipping...' % pdbcode)
            else:
                # Accede to GPCRdb main page of this refined structure
                response = requests.get(refined_baselink+pdbcode)
                if response.ok:
                    # Find the link to the refined structure and donwload it
                    soup = BeautifulSoup(response.content, 'html.parser')
                    btn = soup.find('a',attrs={'id':'download_btn'})
                    btn1 = soup.find('a',attrs={'id':'download_btn1'})
                    link = btn.get('href') if btn else btn1.get('href')
                    response = requests.get(refined_baselink+link)
                    pdbdata = response.content
                    with open(outpdb,'wb') as out:
                        out.write(pdbdata)
                    if response.ok:
                        #Extract pdb file from download zip and save it into a separate file
                        zippy = zipfile.ZipFile(io.BytesIO(response.content))
                        for zippedfile in zippy.namelist():
                            if zippedfile.endswith('.pdb'):
                                pdbfile = zippy.open(zippedfile)
                                content = pdbfile.read()
                                with open(outpdb,'wb') as out:
                                    out.write(content)
                    else:
                        print('could not find link to %s refined structure. Skipping...' % (pdbcode))
                else:
                    print('could not download %s refined structure. Skipping...' % (pdbcode))
            
        except Exception as E:
            print("something failed in downloading refined structure of %s: %s"%(pdbcode,E))

def get_GPCRnames(pdbset):
    """
    Get receptor name from pdbcode via GPCRdb
    """
    GPCRdata = []
    for pdbcode in pdbset:
        pdbdict = requests.get("https://gpcrdb.org/services/structure/"+pdbcode).json()
        GPCRdata.append("%s,%s,%s"%(pdbcode,pdbdict['protein'], pdbdict['signalling_protein']['data']['entity1']['entry_name']))

    return GPCRdata

def create_files(pdbfile, smalmol, name, smalmol_folder, hydrogenate_ligand=True):

    #filenames
    smalmol_pdb = "%s%s_%s.pdb"%(smalmol_folder,name,smalmol) 
    smalmol_mol2 = "%s%s_%s.mol2"%(smalmol_folder,name,smalmol)
    smalmol_script = "%s%s_%s.py"%(smalmol_folder,name,smalmol)
    
    # Skip if mol2file already exists
    if os.path.exists(smalmol_mol2):
        return
        pass

    # Load full structure of system
    mymol = Molecule(pdbfile, validateElements=False)


    # Get resid of first instance of smalmol (we only need one molecule of smaloml, regardless of how many there are in the system)
    lig_resids = mymol.get('resid', 'resname "%s"'%smalmol)
    if not len(lig_resids):
        print("Ligand %s not found in system %s. Skiping"%(smalmol, name))
        return
    resid = lig_resids[0]

    # Write PDB file of smalmol 
    os.makedirs(smalmol_folder, exist_ok=True)
    mymol.write(smalmol_pdb, 'resname %s and resid %s'%(smalmol,str(resid)))

    # Convert it to mol2file, and add Hydrogens to molecule (using chimera)
    if hydrogenate_ligand:
        chimera_addH(smalmol_pdb, smalmol_mol2, smalmol_script, smalmol)
    else:
        chimera_conversion(smalmol_pdb, smalmol_mol2, smalmol_script, smalmol)

def ligand_dictionary(struc_dict, ligandsdict_path, modres_path, basepath, blacklist = {}, hydrogenate_ligands=True):
    """
    Create dictionary with ligand names and ligand ResNames of each of the structures we need to simulate,
    and store the resutls in a json file
    """
    # Read existing ligands dictionary, or create it if it does not exists yet 
    if os.path.exists(ligandsdict_path):
        ligandsdict = json_dict(ligandsdict_path)
    else:
        ligandsdict = {}

    # Read existing modified residues dictionary, or create it if it does not exists yet 
    if os.path.exists(modres_path):
        modresdict = json_dict(modres_path)
    else:
        modresdict = {}
        
    # Iterate over non-yet-simulated structures, and get their ligand information from rcsb (PDB's web api)    
    for (pdb_code,pdbfile) in struc_dict.items():
        #Do not repeat simulations
        if pdb_code in ligandsdict:
            continue
        else:
            ligandsdict[pdb_code] = {}
            # Extract ligand information from rcsb using api 
            reponse_dict = requests.get('https://data.rcsb.org/graphql?query=\
                        {entry(entry_id: "'+pdb_code+'"){\
                            rcsb_entry_info{\
                                nonpolymer_bound_components\
                            }\
                            nonpolymer_entities{\
                                pdbx_entity_nonpoly{comp_id,name}\
                            }\
                            polymer_entities{\
                                entity_poly{rcsb_non_std_monomers}\
                            }\
                        }\
                     }').json()['data']

            #If not an apoform, store obtained ligand information
            if reponse_dict['entry']['nonpolymer_entities']:
                covlig = reponse_dict['entry']['rcsb_entry_info']['nonpolymer_bound_components']
                for ligand in reponse_dict['entry']['nonpolymer_entities']:
                    ligandResname = ligand['pdbx_entity_nonpoly']['comp_id']
                    ligandName = ligand['pdbx_entity_nonpoly']['name']
                    ision = (len(ligandResname) == 2)
                    iscovbound = (ligandResname in covlig) if (covlig) and not (ision) else False
                    # Skip apoforms, blacklisted molecules and ions
                    if (ligandResname == 'null') or (ligandResname in blacklist) or (len(ligandResname)==2):
                        continue
                    else:
                        ligandsdict[pdb_code][ligandResname] = (ligandName, iscovbound)
                        
            # If contains modified residues, store information as well
            for poly in reponse_dict['entry']['polymer_entities']:
                modres = poly['entity_poly']['rcsb_non_std_monomers']
                if modres and len(modres) > 0 :
                    new_modres = []
                    for modre in modres:
                        if modre not in blacklist:
                            new_modres.append(modre)
                    if len(new_modres) > 0:
                        modresdict[pdb_code] = new_modres
                
    # Save dicts as jsons
    with open(ligandsdict_path, 'w') as jsonfile:
        json.dump(ligandsdict, jsonfile, ensure_ascii=False, indent = 4)            
    with open(modres_path, 'w') as jsonfile:
        json.dump(modresdict, jsonfile, ensure_ascii=False, indent = 4)            
    
    # Create ligands set from previou dictionary, and remove blacklisted ligands
    ligandsset = { ligcode  for pdbcode in ligandsdict for ligcode in ligandsdict[pdbcode] }
    ligandsset = ligandsset - blacklist
    
    # Create ligand and modified residues pdb files
    for pdbcode in ligandsdict:
        try:

            # Load curated structure of full system
            pdbfile = struc_dict[pdbcode]
            if not os.path.exists(pdbfile):
                print("Curated file for %s not avaliable. Ligands could not be obtained"%pdbcode)
                continue

            # For each regular ligand, create a mol2 file
            for lig in ligandsdict[pdbcode]:
                if lig not in noparams_ligs: # There are standard parameters for these two
                    is_covlig = ligandsdict[pdbcode][lig][1]
                    # Create mol2 files for regular ligands. Covalent ones go separatedly
                    if not is_covlig: 
                        create_files(pdbfile, lig, pdbcode, '%stoppar/Ligands/%s/'%(basepath,lig),  hydrogenate_ligands)

        except FileNotFoundError as e:
            print("Error: %s"%e)

    return(ligandsdict, ligandsset, modresdict)

def find_aminergic(ligandsdict, gpcrdb_dict, blacklist, strucpath):
    """
    Identify which of the ligands present in this dict are associated to aminergic receptors in GPCRdb
    """
    amiligs = {}
    for pdbcode in ligandsdict:
        # If this pdbcode is a classA (001) aminergic (001)
        if gpcrdb_dict[pdbcode]['family'].startswith('001_001'):
            
            #Load system of pdbcode (if exists)
            pdbfiles = glob(str("%s%s/*.pdb" % (strucpath,pdbcode)))
            if any(pdbfiles):
                mol = Molecule(pdbfiles[0], validateElements=False)
        
            # For each ligand in this system
            for lig in ligandsdict[pdbcode]:
                if (lig not in blacklist) and (len(lig)>2):
                    # Check if there is any N atom in the ligand making a H-bond with an asp/glm
                    natom = mol.get('name', 'resname '+lig+' and element N and (within 3.5 of (resname ASP GLU))')
                    if any(natom):
                        amiligs[lig] = list(natom)
            
    return amiligs

def run_propka(infile, outname, strucpath):
    """
    Run propka on a pdbfile (infile) to stablish ionization of residues
    """

    import propka
    pdbfile = strucpath+outname+'.pdb'
    # Skip if already propkaed
    if os.path.exists(pdbfile):
        return Molecule(pdbfile, validateElements=False)
    
    # Run propka
    mol = propka.run.single(infile)
    #Convert pka to pdb (more or less)
    outfile = open(pdbfile,'w') 
    with open(glob('*.propka_input')[0],'r') as pkafile:
        for line in pkafile:
            outfile.write(line[:55]+'\n')
    outfile.close()
    
    # Delete files created by propka
    for a in glob('*propka_input')+glob('*pka'):
        os.remove(a)

    return Molecule(pdbfile, validateElements=False)

def chimera_addH_pdb(infile, outfile, pyfile, resname="COV"):
    """
    Add Hydrogens to PDB, and save again as PDB
    """
    with open(pyfile,'w') as f:
        f.write("""
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
os.chdir(".") 
rc("open %s") 
rc("addh")
rc("write 0 %s")
""" % (infile, outfile))

    # Execute script in chimera 
    os.system('chimera --nogui %s'%pyfile)


def chimera_conversion(infile, mol2file, pyfile, resname="COV"):
    """
    Write a python scripts with instructions for chimera to read a pdb and save a mol2file, without further changes
    """
    with open(pyfile,'w') as f:
        f.write("""
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open %s") 
rc("setattr m name '%s'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "%s") 
""" % (infile, resname, mol2file))

    # Execute script in chimera 
    os.system('chimera --nogui %s'%pyfile)

def chimera_addH(infile, mol2file, pyfile, resname="COV"):
    """
    # Write a python script giving instructions to chimera for adding hydrogens to our molecule
    # and converting it to mol2 file (for paramchem, you know)
    """
    
    with open(pyfile,'w') as f:
        f.write("""
import os 
from chimera import runCommand as rc 
from chimera import replyobj,openModels,Molecule 
from WriteMol2 import writeMol2 
os.chdir(".") 
rc("open %s") 
rc("addh") 
rc("setattr m name '%s'") 
writeMol2(openModels.list(modelTypes=[Molecule]), "%s") 
""" % (infile, resname, mol2file))

    # Execute script in chimera 
    os.system('chimera --nogui %s'%pyfile)

def save_smalmol_mol2(input_dict, basepath, hydrogenate_ligands=True):
    """
    From the submitted PDB files, extract PDB of small molecule and convert it into mol2
    """
    
    # For everyone of the files we want to simulate
    ligandsdict = {}
    for entry in input_dict:
        pdbfile = basepath+entry['pdbfile']
        ligandsdict[entry['name']] = {}
        # For every small molecule in this system
        for smalmol in entry['ligands']:
            ligandsdict[entry['name']][smalmol['resname']] = [smalmol['name'],smalmol['covalently_bound']]
            create_files(pdbfile,smalmol['resname'],entry['name'], '%stoppar/Ligands/%s/'%(basepath,smalmol['resname']), hydrogenate_ligands)
        # For every modified residue in this system
        for modres in entry['modres']:
            create_files(pdbfile,smalmol,entry['name'], '%stoppar/mod_residues/%s/'%(basepath,modres), hydrogenate_ligands)

    modresdict= { a['name'] : a['modres'] for a in input_dict}
    pdbfilesdict= { a['name'] : a['pdbfile'] for a in input_dict}
    return (modresdict,ligandsdict, pdbfilesdict)

def download_structures(pdbcode, outfile):
    """
    Download PDB structure of certain pdbcode
    """

    reponse = requests.get('https://files.rcsb.org/download/'+pdbcode+'.pdb')
    if reponse.ok:
        with open(outfile,'wb') as infile:
            infile.write(reponse.content)
    else:
        print("No PDB file avalible for "+pdbcode+" in rcsb.org")

"""def paramchem(username,password,ligpath,ligcode,pdbcode):
    # Submit ligand to paramchem to obtain topology-parametrs file OBSOLETE
    
    # Make two paramchems: one with legacy parameters and the other with modern ones (for hallogen systems)
    for c_option in ['c','']:
        
        #Pattern to find non-HTMD-compatible lines
        lph_pat = re.compile('^ATOM.*LPH|LONEPAIR')

        # Ligand mol2file
        mol2file = ligpath+pdbcode+'_'+ligcode+".mol2"
        # topology-parametetrs  and warnings file output
        if c_option:
            topparfile_path ="%s%s_legacy_toppar.str"%(ligpath,pdbcode)
            erwar_path ="%s%s_legacy_paramchem_stder.txt"%(ligpath,pdbcode)
        else:
            topparfile_path ="%s%s_toppar.str"%(ligpath,pdbcode)
            erwar_path ="%s%s_paramchem_stder.txt"%(ligpath,pdbcode)
            
        # Skip if input mol2 file does not exists
        if not os.path.exists(mol2file):        
            print("file %s does not exist. Skipping...."%mol2file)
            continue
        # Omit ligand if its toppar file already exists
        if os.path.exists(topparfile_path):
            leg = 'legacy' if c_option else 'latest'
            print('%s toppar for ligand %s already exists. Skipping...' % (leg,ligcode))
            continue
        
        # Ligands for which chimera creates wrong topology have to be autodetermined in paramchem
        badtop_lig = {'7AB'}
        b_option = 'b' if ligcode in badtop_lig else ''

        #Open ligandfile to upload in paramchem
        myfile = {
                'filename': open(mol2file)
        }

        # Define POST variables for login in Paramchem
        datalogin = {
            'usrName': username,
            'curPwd': password,
            'rndNo': str(random.randrange(100000000, 999999999)),
            'submitBtn': 'Submit',
        }

        # Variables for paramchem
        myparam = {
                #'param_a': 'a', #Include parameters usually included in newer versions of CGenff (versions that we cant use)
                'param_b' : b_option , # Make paramchem determine bond order (active if covalent-bound ligand)
                'param_c': c_option, # Use CGenFF legacy v1.0, for HTMD is not yet prepared for newer Charmm versions             
        }

        # Open web session
        with requests.Session() as s:

            # Login into paramchem
            response_login = s.post('https://cgenff.silcsbio.com/userAccount/userWelcome.php', 
                       data=datalogin,
                       verify=False)

            # Submit our ligand molecule into paramchem
            response_upload = s.post('https://cgenff.silcsbio.com/initguess/processdata.php', 
                       files = myfile,
                       data = myparam,
                        )

            # Download Topology-parameters of our molecule file from paramchem, and store it.
            # But remember submissions in paramchem are limited weekly
            # Download also stderr, just in case
            match = re.search('<path>(\w+)</path>', response_upload.text)
            if match:
                code = match.group(1)
                response_ligfile = s.get('https://cgenff.silcsbio.com/initguess/filedownload.php?file=%s/%s_%s.str'%(code,pdbcode,ligcode))
                response_stder = s.get('https://cgenff.silcsbio.com/initguess/filedownload.php?file=%s/%s_%s.err'%(code,pdbcode,ligcode))
                with open(topparfile_path, 'wb') as topparfile:
                        topparfile.write(response_ligfile.content)
                with open(erwar_path, 'wb') as erwar:
                        erwar.write(response_stder.content)                            
            else:
                print('Your paramchem account has reached its weekly submission limit. Please, intrudce a new account or wait to the next monday to continue')            
                return

            #Delete lines with LPH (new feature from CHARMM not tolerated by HTMD)
            with open(topparfile_path, "r") as f:
                lines = f.readlines()
            with open(topparfile_path, "w") as f:
                for line in lines:
                    if not re.match(lph_pat, line):
                        f.write(line)
"""        

def get_params_cgenff(username,password,mol2path,outfolder):
    """
    Obtain parameters for molecule in mol2path using CGenFF app web service
    """

    # URLs for the Identity Toolkit API
    API_KEY = 'AIzaSyA_k5QDKh2dAZLzcg9MsG4tB_BfSIiEGJw'
    SIGN_IN_URL = 'https://identitytoolkit.googleapis.com/v1/accounts:signInWithPassword?key='+API_KEY
    SECURE_RESOURCE_URL = 'https://app.cgenff.com/homepage'  # Replace with the actual secure resource URL

    def sign_in(email, password):
        payload = {
            'email': email,
            'password': password,
            'returnSecureToken': True
        }
        response = requests.post(SIGN_IN_URL, json=payload)
        response.raise_for_status()  # Raise an error if the request was unsuccessful
        return response.json()

    def submit_molecule(file_path, session, user_id, job_id):
        """Submit mol2 file into CgenFF app"""

        data = {
            'userId' : user_id,
            'jobIDFolder' : job_id
        }

        headers = {
            'Authorization': f'Bearer {id_token}'
        }

        # Upload the molecule file
        files = {'mol2File': open(file_path, 'rb')}
        response = session.post('https://app.cgenff.com/jobs/file-management/uploadMol2File', files=files, headers=headers, data=data)
        response.raise_for_status()
        
        # Return the task ID or results URL if provided
        print('submited_mol',response)

    def submit_cgenff(session, user_id, job_id, stem):
        """
        Process molecule in CGenFF to obtain parameters
        """

        print('submitting into cgenff....')
        headers = {
            'User-Agent': 'python-requests/2.22.0', 
            'Accept-Encoding': 'gzip, deflate',
            'Accept': '*/*', 
            'Connection': 'keep-alive', 
            'Content-Length': '2705',
            'Content-Type': 'application/json'
            }

        data = json.dumps({
            "userId":user_id,
            "jobId":job_id,
            "stem":stem,
            "options":{"mol2out":True,"warning":False,"debug":False,"copyParams":False}
            })
        response = session.post('https://app.cgenff.com/job-managers/submitCGenFF', headers=headers, data=data)
        response.raise_for_status()

    def str_cgenff(session, user_id, job_id, stem, outfolder='./'):
        """
        Download STR file with topologies and parameters
        """

        data = {
            "userId":user_id,
            "jobId":job_id,
            "stem":stem,
        }

        headers = {
            'User-Agent': 'python-requests/2.22.0', 
            'Accept-Encoding': 'gzip, deflate',
            'Accept': '*/*', 
            'Connection': 'keep-alive', 
            'Content-Length': '2705',
            }

        # Download both error log file (.err) and the actual .str file with parameters and topology of molecule
        for filext in ('.err','.str'):
            response = session.get('https://app.cgenff.com/build/%s/%s/%s%s'%(data['userId'],data['jobId'],stem,filext), headers=headers, data=data)
            response.raise_for_status()

            # Write returned parameters in a file
            outfile = outfolder+stem+'_toppar'+filext
            with open(outfile, 'wb') as filehand:
                filehand.write(response.content)
        return(outfile)

    # Extrac filename without folders or filextension (requried for CGenFF parameters)
    stem = os.path.splitext(os.path.basename(mol2path))[0]

    # Skip if molecule has no mol2
    if not os.path.exists(mol2path):
        print("file %s does not exist. Skipping...."%mol2path)
        return

    # Create a session to manage cookies
    session = requests.Session()

    # Sign in to get the ID token
    auth_response = sign_in(username, password)
    id_token = auth_response['idToken']
    user_id = auth_response['localId']

    # Create a random job id
    job_id = str(random.randrange(1000000000000, 5555555555555))

    # Submit molecule into CGenFF
    submit_molecule(mol2path,session, user_id, job_id)

    # Proces molecule parameters
    submit_cgenff(session, user_id, job_id,stem)

    # Download said parameters
    outfile = str_cgenff(session, user_id, job_id, stem, outfolder)

    #Delete lines with LPH (feature from CHARMM not tolerated by HTMD)
    # lph_pat = re.compile('^ATOM.*LPH|LONEPAIR')
    # with open(outfile, "r") as f:
    #     lines = f.readlines()
    # with open(outfile, "w") as f:
    #     for line in lines:
    #         if not re.match(lph_pat, line):
    #             f.write(line)

def get_lig_toppar(ligandsdict, basepath, username, password, pdbfiles = {}):
    """
    Get the topology-parameters string file from paramchem for the submited ligand PDB codes
    ALERT: paramchem only allows 100 submissions by month, so it may be possible that not all 
    parameters are obtained
    """
    
    # Iterate over ligands
    for pdbcode in ligandsdict:
        for ligcode in ligandsdict[pdbcode]:
            
            # Skip Retinol and cholesterol (we already have parameters for these ones) or ions (unprocessable for paramchem) 
            # or blacklisted molecules
            if (ligcode in noparams_ligs) or (len(ligcode)<3) or (ligcode in detergent_blacklist.union(glucids_blacklist)):
                continue
            
            # Stablish folder for ligand parameters (depending on covalent-bound ligand or not)
            ligcov = ligandsdict[pdbcode][ligcode][1]
            ligpath = basepath+"toppar/Ligands/"+ligcode+'/'
            mol2file = ligpath+pdbcode+'_'+ligcode+'.mol2'
            strfile = "%s%s_%s_toppar.str"%(ligpath,pdbcode,ligcode)

            # Skip already-parameterized ligands
            if os.path.exists(strfile):
                continue
            
            print('Getting toppar file for ligand %s ' % (ligcode))
                                
            # IF is a covalent-bound ligand
            if ligcov:
                modres_covlig_toppar(pdbcode, ligcode, ligpath, username, password, covlig = True, pdbfiles = pdbfiles)
            else:
                # paramchem(username,password,ligpath,ligcode,pdbcode)
                get_params_cgenff(username,password,mol2file,ligpath)

def fix_ARG(inpath, outpath):
    """
    Fix arginine's (and similars) topology. Chimera makes CZ to have three double bonds, which
    result in crashing paramchem (hypervalent carbon).
    Here, I force any CZ to have single-bonds with NH2 and NH1
    Also, I fix similar problem with -COO_ groups 
    """
    infile = open(inpath, 'r')
    outfile = open(outpath, 'w')
    atoms = False
    bonds = False
    cz = set()
    ne = set()
    c = set()
    oxt = set()
    for line in infile:
        if line.startswith('@<TRIPOS>ATOM'):
            atoms = True
            outfile.write(line)
            continue            
        elif line.startswith('@<TRIPOS>BOND'):
            bonds = True
            atoms = False
            outfile.write(line)
            continue
        elif line.startswith('@<TRIPOS>SUBSTRUCTURE'):
            bonds = False
        if atoms:
            line_list = line.split()
            if line_list[1] == 'CZ':
                cz.add(line_list[0])
                line = line.replace('C.cat', 'C.2  ')
            elif (line_list[1] == 'NH2') or (line_list[1] == 'NE'):
                ne.add(line_list[0])
                line = line.replace('N.3  ', 'N.pl3')
            elif line_list[1] == 'OXT':
                oxt.add(line_list[0])
            elif line_list[1] == 'CYX':
                c.add(line_list[0])
        elif bonds:
            line_list = line.split()
            if (line_list[1] in cz and line_list[2] in ne) or\
            (line_list[1] in ne and line_list[2] in cz) or \
            (line_list[1] in c and line_list[2] in oxt) or \
            (line_list[1] in oxt and line_list[2] in c):
                line_list[3] = '1'
                line = "%s%s%s %s\n"%(line_list[0].rjust(6),line_list[1].rjust(5),
                                  line_list[2].rjust(5),line_list[3])
        outfile.write(line)
    outfile.close()
    infile.close()
    
def get_modres_toppar(modresdict, basepath, username, password, pdbfiles={}):
    """
    Get topologies and parameters for modified residues
    """
    
    for pdbcode in modresdict:
    # for pdbcode in ['5VBL']:

        # For each PDB structure, iterate on its modified residues
        for modres in modresdict[pdbcode]:
            modrespath = basepath+'toppar/mod_residues/'+modres+'/'
            modres_covlig_toppar(pdbcode, modres, modrespath, username, password, pdbfiles)


def renumber_covlig(mol, lig):
    """
    # Set new atom names for ligand (to avoid repeated names with residue)
    """

    for ele in {'N','C','O','S'}:
        i=0
        for atom in mol.get('index',str('resname %s and element %s' %(lig, ele))):
            mol.set('name', ele+str(i), str('index %i and element %s' %(atom, ele)))
            i+=1

    return(mol)


def modres_covlig_toppar(pdbcode, modres, modrespath, username, password, covlig = False, pdbfiles = {}):
    """
    Prepare residue to be parameterized, send it to paramchem, and adjust its topology
    """

    # FIlenames and directories
    mol2chim = modrespath+modres+'_chim.mol2'
    mol2file = modrespath+pdbcode+'_'+modres+'.mol2'
    topparfile = modrespath+pdbcode+'_'+modres+'_toppar.str'
    pdbfile = modrespath+modres+'.pdb'
    pyfile = modrespath+modres+'.py'
    os.makedirs(modrespath, exist_ok = True)

    # Skip existing ones
    if os.path.exists(topparfile):
        print("toppar file for %s in %s already exists. Skipping"%(modres, pdbcode))
        return
        pass

    # Load main structure
    if (pdbcode not in pdbfiles):
        print('pdbfile for %s not avaliable. SKipping...'%pdbcode)
        return
    elif not os.path.exists(pdbfiles[pdbcode]):
        print('mol2file for %s not avaliable. SKipping...'%pdbfiles[pdbcode])
        return
    main_pdb = Molecule(pdbfiles[pdbcode], validateElements=False)

    # Find out if this structure has hallogens
    has_halo = bool(len(main_pdb.get('name', 'element Cl Br I and not ion')))
    
    # To avoid parametrization problems, rename repeated molecules of our modres (aka: only one modres molecule per system)
    modres_resid = main_pdb.get('resid', 'resname '+modres)[0]
    main_pdb.set('resname', 'XYX', 'resname %s and not resid %d'%(modres, modres_resid))

    # If covalent ligand, change resname of ligand+boundresidue to COV
    if covlig:
        main_pdb = renumber_covlig(main_pdb, modres)
        sel = 'resname %s or same residue as (protein and within 2 of (resname %s))' % (modres, modres)
        new_resid = main_pdb.get('resid', 'protein and within 2 of (resname %s)'%modres)[0]
        new_chain = main_pdb.get('chain', 'protein and within 2 of (resname %s)'%modres)[0]
        main_pdb.set("resname", "COV", sel)
        main_pdb.set("resid", str(new_resid), "resname COV")
        main_pdb.set("chain", str(new_chain), "resname COV")
        main_pdb.set("segid", 'COV', "resname COV")
        modres = "COV"
    
    # Save in a list the residue's binding atoms
    chains = main_pdb.get('chain', 'resname '+modres)
    if len(chains) > 0:
        modres_chain = chains[0]
    else:
        print('modified residue '+modres+' not present in structure '+pdbcode+'. Skipping...')
        return
    

    # In the main structure, delete everything except for the ligand, neighboring residues
    # and the atoms bounding these neighboring residues
    modres_resid = main_pdb.get('resid', 'resname '+modres)[0]
    preserve_resids_list = [modres_resid, modres_resid-1, modres_resid+1, modres_resid-2, modres_resid+2]
    preserve_resids_list = [item for item in preserve_resids_list if item >= 0]
    preserve_resids =' '.join(map(str, preserve_resids_list))
    main_pdb.filter('(resid %s) and chain %s'% (preserve_resids, modres_chain))

    # Iterate over modres atoms to find the atoms from neighboring residues bound to modres
    boundatoms = ''
    modres_boundatoms = set() 
    replace_atoms = dict()
    i=0
    for atom in main_pdb.get('index','resname '+modres):
        atom_name = main_pdb.get('name', 'index '+str(atom))[0]
        for nearby_atom in main_pdb.get('index', 'not resname %s and within 2 of (index %s)'%(modres,atom)):
            modres_boundatoms.add('ATOM '+atom_name+' ')
            n_atom_name = main_pdb.get('name', 'index '+str(nearby_atom))[0]
            # Skip bonds modresN -> nonmodresC (Repeated bonds give problems in PSF files)
            if (n_atom_name.startswith('C') and (n_atom_name != 'C')) or \
            ((n_atom_name == 'C') and (atom_name == 'N')):
                continue
            else:
                yname = 'Y%d'%i
                main_pdb.set('name', yname, 'index '+str(nearby_atom))
                i += 1
            if n_atom_name.startswith('C'):
                boundatoms += yname+' '
                replace_atoms[' '+yname+' '] = ' -'+n_atom_name+' '
            if n_atom_name.startswith('N'):
                boundatoms += yname+' '
                replace_atoms[' '+yname+' '] = ' +'+n_atom_name+' '

    # 'C'-named atoms are not properly hidgogenated in chimera
    # Same happens with O-named atoms
    # And ASP and GLU are left unprotonated, which crashes paramchem
    main_pdb.set('name','CYX', 'name C')
    main_pdb.set('name','OYX', 'name O')
    main_pdb.set('resname','GLH', 'resname GLU')
    main_pdb.set('resname','ASH', 'resname ASP')

    # Set a common segid and write a PDB
    main_pdb.set('segid', 'X')
    main_pdb.write(pdbfile)

    # Add hydrogens with chimera (THIS GIVES A LOT OF PROBLEMS. AVOID USE)
    chimera_addH(pdbfile, mol2chim, pyfile, modres)

    # Fix arginine bug
    fix_ARG(mol2chim,mol2file)
    os.remove(mol2chim)

    # Open the hydrogenated molecule and prepare it for paramchem
    mol2 = Molecule(mol2file, validateElements=False)

    # Set atom names of non-modres (or bound to modres) as XXX
    if boundatoms:
        mol2.set('name', 'X', 'not (resname %s or name %s)'%(modres, boundatoms))
    else:
        mol2.set('name', 'X', 'not resname '+modres)
    mol2.set('name','C', 'name CYX')# Reverse past changes
    mol2.set('name','O', 'name OYX')# Reverse past changes
    mol2.set('resname', modres, 'all')
    mol2.write(mol2file)

    # Use paramchem to obtain toppar file for our ligand-residue thing
    get_params_cgenff(username,password,mol2file,modrespath)
    # paramchem(username,password,modrespath,modres,pdbcode)

    #Rename output files from paramchem
    pre_topparfile = modrespath+'pre_'+pdbcode+'_'+modres+'.str'
    os.rename(topparfile, pre_topparfile)

    # Create new file for parameters
    with open(pre_topparfile,'r') as f:
        topparlines = f.readlines()
    new_f = open(topparfile,'w')

    Xpat = re.compile(" X\w* ")
    typepat = re.compile("ATOM \w+\s+(\w+)")
    for line in topparlines:
        if line.startswith('ATOM'):
            # Skip atom lines refering to atoms from neighbor residues
            if line.startswith(('ATOM X','ATOM Y')):
                continue

        elif line.startswith(('BOND','IMPR')):
            if Xpat.search(line):
                continue
            elif any(a in line for a in replace_atoms):
                for atom in replace_atoms:
                    line = line.replace(atom, replace_atoms[atom])

        new_f.write(line)

    new_f.close()

def hetatm_nucleotides(pdbpath, name):
    """
    Convert all nucleotide residue lines from ATOM into HETATM (they will be excluded 
    by homolwat otherwise).
    """
    mol = Molecule(pdbpath, validateElements=False)
    mol.set('record', 'HETATM', sel='resname ADN GDP GTP')
    mol.write(pdbpath)
    file = open(pdbpath, 'rb')
    return file    

def add_peplig(filename, pdbcode, gpcrdb_dict):
    """
    GPCRdb refined structures of peptide-ligand receptors do not have its corresponding ligand crystalzyed. 
    Therefore, I have to re-add it by superimposing the original PDB structure
    """

    # Skip if not actually a peptide-complex
    if not any([ True for lig in gpcrdb_dict[pdbcode]['ligands'] if lig['type'] in {'protein', 'peptide'} ]):
        return

    # Load my input structure
    mymol = Molecule(filename, validateElements=False)
    # Check if already has PEP ligand for wathever reason
    if len(mymol.get('segid', 'segid PEP')):
        return

    # Load PDB standard structure
    pdbmol = load_mol_pdbcode(pdbcode)
    
    # Align in-preparation molecule system with its original PDB counterpart
    alignment_results = sequenceStructureAlignment(pdbmol, mymol, maxalignments = 1)
    new_pdbmol = alignment_results[0]
    
    # Extract chain names and length from PDB
    chains_dict = requests.get('https://data.rcsb.org/graphql?query=\
        {entry(entry_id: "'+pdbcode+'") {\
            polymer_entities {\
              entity_poly{rcsb_sample_sequence_length, pdbx_strand_id}\
              rcsb_polymer_entity{pdbx_description}\
            }\
          }\
        }\
    ').json()

    # Check which one of this chains corresponds to the ligand protein chain
    ligchain = None
    min_seqlen = 100000000000
    for poly in chains_dict['data']['entry']['polymer_entities']:
        uniname = poly['rcsb_polymer_entity']['pdbx_description']
        uniname = uniname.lower() if uniname else False
        seqlen = int(poly['entity_poly']['rcsb_sample_sequence_length'])
        chainId = poly['entity_poly']['pdbx_strand_id'].split(',')[0]# In case there are twin systems like 4bUO

        # If uniprot names of the chain exist and contain 'receptor' or 'G-protein', then for sure this is is not the ligand chain
        if uniname and any(word in uniname for word in ['receptor','Guanine nucleotide-binding protein']):
            continue

        # Compare lengths and take shortest
        if seqlen < min_seqlen:
            min_seqlen = seqlen
            ligchain = chainId

    # Take ligand chain from the original pdb
    new_pdbmol.filter("chain "+ligchain+" and not water")

    # Change chain id and seg id of peptide ligand
    new_pdbmol.set('segid', 'PEP', ' chain '+ligchain)
    new_pdbmol.set('chain', 'Q', 'chain '+ligchain)

    # Append new ligand to structure
    mymol.append(new_pdbmol)
    
    # Re-write result
    mymol.write(filename)

def needs_sodium(gpcrdb_dict, pdbcode):
    """
    Determine (based on information in GPCRdb) if structures require 2x50 sodium or not
    """
    agoset = {'Agonist', 'Agonist (partial)','Ago-PAM', 'AgonistNAM', 'AgonistPAM', 'Allosteric agonist', 'PAM', 'NAM', 'unknown'}
    ligands = gpcrdb_dict[pdbcode]['ligands']
    family = gpcrdb_dict[pdbcode]['family'] 
    # Gess if this structure needs or not a sodium
    sod=True
    if not (family.startswith('001')) or (family.startswith('001_002_023')) : # If not a class-A gpcr OR is a orexin
        sod = False
    elif (len(ligands) == 0) or (ligands[0]["name"]=="Apo (no ligand)"): # If apoform, do not add the ion (we want to see how it enters)
        sod = False
    elif (ligands[0]['function']  in agoset): # If agonist-bound complex
        sod = False

    return sod


def internal_waters(pdbpath, pdbcode, gpcrdb_dict, apo=False, sod=False,chain=False):
    """
    Place internal waters and E2x50 sodium in GPCR structure using homolwat online tool
    """

    # Set names 
    name = os.path.splitext(pdbpath)[0]

    # Check if there is already a watered structure here
    watered_filename = name+'_apoHW.pdb' if apo else name+'_HW.pdb'
    if os.path.exists(watered_filename) > 0:
        print("Structure %s already has a watered version. Skipping..." % pdbcode)
        return ( sod, watered_filename)
    else:
        print("Adding internal waters to structure")

    #Load pdb file
    pdbfiles = {"file" : hetatm_nucleotides(pdbpath, name)} 

    # Open web session
    with requests.Session() as s:
        # Load our desired structure into homolwat, as a PDB file
        response_loadfile = s.post("https://alf06.uab.es/homolwat/run_HW",
                                   files = pdbfiles
        )
        # Check if we are in the multiple-chain menu
        soup = BeautifulSoup(response_loadfile.text,'html')
        form_action = soup.find('form').get('action')
        
        if form_action=="run_hw_multi":
            # Select GPCR chain (first one, if none has been specified)
            first_chain = soup.find('option').get('value')
            chainid = chain if chain else first_chain
            filename = soup.find('input', attrs={'name':'filename'}).get('value')
            response_loadfile = s.post(
                "https://alf06.uab.es/homolwat/run_hw_multi",
                  data= {
                    'filename' : filename,
                    'option_gpcr' : chainid
                  }
            )
            
        # Get name and number of our query. Required for following steps
        soup = BeautifulSoup(response_loadfile.text,'html')
        query_num = soup.find('input',attrs={'name' : 'query_num'}).get('value')
        query_name = soup.find('input',attrs={'name' : 'query_name'}).get('value')

        # Analyze structure and place internal waters
        response_solvate = s.post("https://alf06.uab.es/homolwat/solvate",
                  data = {
                    "query_name": query_name,
                    "query_num": query_num,
                    "ch_rest": '[]',
                    "option_state": "inac",
                    "option_sodium": "sod_yes" if sod else "sod_no",
                    "option_dowser": "dow_no",
                    "p_ident": ""
                  }
        )

        # Determine name for output file
        query_num_noapo = query_num.split("'")[1]
        if apo:
            apofix = '_apo'
        else:
            apofix = '_'



        # Download results in zip            
        query_name_nofilext = os.path.splitext(query_name)[0]
        response_download = s.post("https://alf06.uab.es/homolwat/download_file/",
              data = {
                    "query_num": query_num_noapo,
                    "query_name": query_name_nofilext
              })

        # Unzip and extract watered strucutre file
        req = response_download.request
        water_lines = []
        if response_download.ok:
            zippy = zipfile.ZipFile(io.BytesIO(response_download.content))
            zippy.extract(query_name_nofilext+"_HW.pdb", path='./')
            # This program sodiums are wrongly formated, so the name of the homolog
            # simulation occupies the charge column by error. Here i solve this
            in_hw = open('./' + query_name_nofilext+"_HW.pdb", 'r')
            out_hw = open(watered_filename, 'w')
            for line in in_hw:
                if (len(line) > 79) and (line[79] != " "):
                    linelist = list(line)
                    linelist[79] = " "
                    line = "".join(linelist)

                # Save water and sodium lines added by Homolwat
                if (len(line)>81) and (line[17:20] in ['HOH',' NA']):
                    water_lines.append(line)

            # Merge homolwat lines to main refined input structure
            with open(pdbpath, 'r') as ref:
                for line in ref:
                    if not line.startswith('END'): 
                        out_hw.write(line)
                for line in water_lines:
                    out_hw.write(line)
                out_hw.write('END')

            os.remove('./' + query_name_nofilext+"_HW.pdb")
            in_hw.close()
            out_hw.close()

        else:
            print("could not add internal waters to %s. Skipping..." % (pdbpath))
            
    return (sod, watered_filename)

def dowser_waters(pdbpath, dowserbin, outpath):
    """
    Use dowser (old version, not "dowser++") to add internal waters to the system
    """
    os.system(dowserbin+" "+pdbpath)
    
    # Merge obtained waters with system pdb
    pdbfile = open(pdbpath,'r')
    waterfile = open('dowserwat_all.pdb','r')
    conectors = []
    with open(outpath,'w') as outfile:
        for file in [pdbfile, waterfile]:
            for line in file:
                # These must go at the end
                if line.startswith('CONECT'):
                    conectors.append(line)
                # This always give problems. Skip them
                elif not any(a in line for a in ['END','TER']):
                    outfile.write(line)

        # Append topology lines            
        for line in conectors:
            outfile.write(line)

    pdbfile.close()
    waterfile.close()
    
    # Remove remaining dowser files
    for file in ['dowserwat_all.pdb', 'dowserwat.pdb', 'intsurf.pdb','reform.pdb']:
        os.remove(file)
        

def remove_ligmols(gpcrdb_mol, blacklist, gpcrdb_dict, pdbcode, apo):
    """
    Check some stuff concerning to ligand molecules
    """
    # Remove unnecessary ligand molecules: mostly crystalization detergents, quelants, buffers,
    # or post-traductional glicosilations
    gpcrdb_mol.remove('resname '+' '.join(blacklist))

    # Remove 2x50Sodium from non-A-class GPCRs
    if not gpcrdb_dict[pdbcode]['family'].startswith('001'):
        gpcrdb_mol.remove('element NA')

    # Remove ligands if is apo form
    if apo:
        gpcrdb_mol.remove('not ((same chain as protein) or lipid or ion)')

    return gpcrdb_mol

def get_thickness(pdbcode):
    """
    Return membrenae thickness for this pdbcode according to OPM database
    """
    
    # Search pdb code in OPM database, and take thickness and opm_id
    searchlink = 'https://opm-back.cc.lehigh.edu/opm-backend//primary_structures?search='+pdbcode+'&sort=&pageSize=100'
    response_dict=eval(requests.get(searchlink).content.decode('UTF-8')\
                       .replace('true', 'True')\
                       .replace('false','False')\
                       .replace('null','None'))

    # If no thickness avalible for this system, get the one of a random GPCR (not that much of a difference anyway)
    if len(response_dict['objects']):
        thickness = response_dict['objects'][0]['thickness']
    else: 
        thickness = get_thickness('4EJ4')

    return(thickness)

def get_opm(pdbcode):
    """
    Download and load OPM molecule for this PDB code. 
    If there is no opm structure for this pdbcode, download the substitute structure OPM seems fit
    Return thickness and opm molecule
    """
    
    # Search pdb code in OPM database, and take thickness and opm_id
    searchlink = 'https://opm-back.cc.lehigh.edu/opm-backend//primary_structures?search='+pdbcode+'&sort=&pageSize=100'
    headers = {
        'Accept' : 'application/json, text/plain, */*', 
        'Accept-Encoding' : 'gzip, deflate, br, zstd', 
        'Origin' : 'https://opm.phar.umich.edu', 
        'Referer' : 'https://opm.phar.umich.edu/'
        }
    response_dict=requests.get(searchlink, headers=headers).json()

    # If the PDBcode has no avaliable OPM structure, use a random pdbcode as reference (4EJ4 in this case)
    if not len(response_dict['objects']):
        searchlink = 'https://opm-back.cc.lehigh.edu/opm-backend//primary_structures?search=4EJ4&sort=&pageSize=100'
        response_dict=requests.get(searchlink, headers=headers).json()
    
    thickness = response_dict['objects'][0]['thickness']
    opm_id = response_dict['objects'][0]['id']

    # Throught opm id, get PDB id of the reference structure in OPM for this pdbcode
    opmlink = 'https://opm-back.cc.lehigh.edu/opm-backend//primary_structures/'+str(opm_id)
    response_dict=requests.get(opmlink, headers=headers).json()
    new_pdbcode = response_dict['pdbid'].lower()
    
    # Download OPM pdb file, and save all lines not involving a DUM residue in a temp file 
    response_pdb = requests.get('https://opm-assets.storage.googleapis.com/pdb/'+new_pdbcode+'.pdb', headers=headers)
    line_list = (response_pdb.content.decode('UTF-8').split('\n'))
    tmpout = tempfile.NamedTemporaryFile(mode='w',suffix='.pdb')
    for line in line_list:
        if ' DUM ' not in line:
            tmpout.write(line+'\n')
    name = tmpout.name
    mol = Molecule(name, validateElements=False)
    tmpout.close()
    
    return (thickness, mol)

def remove_artifacts(pdbcode, mol, ligdict, accepted_ligdict):
    """
    Remove any small molecules included in ligdict but not in accepted_ligdict.
    The intention is to remove unnecessary small molecules, like detergents 
    or ligands from removed parts of the protein
    """
    tofilter = ""
    for smalmol in ligdict[pdbcode]:
        if smalmol not in accepted_ligdict[pdbcode]:
            print('no '+smalmol)
            tofilter += smalmol + " "
    if tofilter:
        mol.filter('not resname '+tofilter)
    return mol
                            
#####################
## Ismael's Functions
#####################

def renumber_resid_vmd(mol,sel,by=3,start=1):
    tmpin = tempfile.NamedTemporaryFile(suffix='.pdb')
    mol.write(tmpin.name)
    viewer = getCurrentViewer(dispdev='text')
    viewer.send('set molid [mol new {%s}]' % tmpin.name)
    tmpin.close()
    tmpout = tempfile.NamedTemporaryFile(suffix='.pdb')
    
    # What a wierd way of deciding if "by_segid", "by_resname" or by both
    option_num = 2
    max_value = 2**option_num - 1
    if by > max_value:
        raise ValueError('Maximum value for "by" keyword is %d.' % max_value)
    if by < 1:
        raise ValueError('Minimum value for "by" keyword is "1".')
    bin_by = format(by,'0'+str(option_num)+'b')
    option_array = [bool(int(i)) for i in bin_by]
    by_segid = option_array[0]
    by_resname = option_array[1]
    
    if by_segid:
        segids = set(mol.get('segid',sel=sel))      
        for segid in segids:
            if by_resname:
                resnames = set(mol.get('resname',sel='(%s) and segid %s' % (sel,segid)))
                for resname in resnames:
                    lsel = '(%s) and (segid %s) and (resname %s)' % (sel,segid,resname)
                    viewer = renumber_resid_by_resid_vmd(lsel,mol,viewer,start=start)
            else:
                lsel = '(%s) and (segid %s)' % (sel,segid)
                viewer = renumber_resid_by_resid_vmd(lsel,mol,viewer,start=start)
    else:                      
        resnames = set(mol.get('resname',sel=sel))
        for resname in resnames:
            lsel = '(%s) and (resname %s)' % (sel,resname)
            viewer = renumber_resid_by_resid_vmd(lsel,mol,viewer,start=start)
    viewer.send('animate write pdb {%s} waitfor all top;exit' % tmpout.name)
    newmol = Molecule(tmpout.name, validateElements=False)
    tmpout.close()
    return newmol

def renumber_resid_by_resid_vmd(sel,mol,viewer,start=1):
    resids = sorted(list(set(mol.get('resid',sel=sel))))
    resids = [str(i) for i in resids]
    viewer.send('proc renum_resid {molid} {set newresid %d; set resids {%s};' % (start,' '.join(resids)) + \
                'set asall [atomselect $molid [concat {(%s) and resid } $resids]];' % sel + \
                '$asall set user 1.00;' + \
                'foreach resid $resids {' + \
                'set as [atomselect $molid [concat {user 1.00 and (%s) and resid } $resid]];' % sel + \
                '$as set resid $newresid; $as set user 0.00; incr newresid}};'+'renum_resid $molid')
    return viewer

def ordered_unique(seq):
    seen = set()
    return [x for x in seq if not (x in seen or seen.add(x))]

def renumber_segments(inputmol,segids,prefix):
    sel = 'segid '+' '.join(segids)
    segids = ordered_unique(inputmol.get('segid',sel=sel))
    if prefix in segids:
        raise ValueError('segid %s already exists.' % prefix)
    
    mol = renumber_resid_vmd(inputmol,sel,by=2)
    # change first segid segment as it is properly renumbered already
    mol.set('segid',prefix,sel='segid '+segids[0])

    if len(segids) > 1:
        # initialize variables for second segid
        curr_segid = prefix
        # get last resid for first segid
        idx_curr_segid = mol.atomselect('segid '+curr_segid)
        prev_resid = len(set(mol.resid[idx_curr_segid]))
        k = 0
        for segid in segids[1:]:
            
            # get last current resid
            idx_segid = mol.atomselect('segid '+segid)
            curr_resid = len(set(mol.resid[idx_segid])) + prev_resid
            if curr_resid <= 9999:
                # join segments resuming the previous resid numbering
                mol = renumber_resid_vmd(mol,'segid '+segid,start=prev_resid+1,by=2)
                mol.segid[idx_segid] = curr_segid
                # get last resid of the current segid for the next loop iteration
                prev_resid = curr_resid
            else:

                # join segments resuming the previous resid numbering up to resid 9999
                sel1 = 'segid '+segid+' and resid <= '+str(9999-prev_resid)
                mol = renumber_resid_vmd(mol,sel1,start=prev_resid+1,by=2)
                mol.set('segid',curr_segid,sel=sel1)
                # define next new segment with resids > 9999
                k +=1
                curr_segid = prefix+str(k)
                if curr_segid in segids:
                    raise ValueError('segid %s already exists.' % curr_segid)
                # resid <= 9999 still preserve the old segid
                idx_curr_segid = mol.atomselect('segid '+segid)
                mol.segid[idx_curr_segid] = curr_segid
                # get last resid of the current segid for the next loop iteration
                mol = renumber_resid_vmd(mol,'segid '+curr_segid,by=2)
                prev_resid = len(set(mol.resid[idx_curr_segid]))
            
        if k > 0:
            if prefix+str(0) in segids:
                print('WARNING: segid %s already exists, using %s instead.' % (prefix,prefix+str(0)))
            else:
                mol.segid[mol.segid == prefix] = prefix+str(0)
        
    return mol

def renumber_resid_by_resid(sel,mol,ordered=False):
    resids = list(set(mol.get('resid',sel=sel)))
    if ordered:
        resids = sorted(resids)
    newresid = 1
    for resid in resids:
        mol.set('resid',newresid,sel='(%s) and (resid %s)' % (sel,resid))
        newresid += 1
    return mol

def renumber_resid(mol,sel,by=3):
    
    # Long story short: by=1 -> by_resname; by=2 -> by_segid; by=3 -> by_segid and by_resname
    # WTF!!!!
    option_num = 2
    max_value = 2**option_num - 1
    if by > max_value:
        raise ValueError('Maximum value for "by" keyword is %d.' % max_value)
    if by < 1:
        raise ValueError('Minimum value for "by" keyword is "1".')
    bin_by = format(by,'0'+str(option_num)+'b')
    option_array = [bool(int(i)) for i in bin_by]
    by_segid = option_array[0]
    by_resname = option_array[1]
    
    if by_segid:
        segids = set(mol.get('segid',sel=sel))      
        for segid in segids:
            if by_resname:
                resnames = set(mol.get('resname',sel='(%s) and segid %s' % (sel,segid)))
                for resname in resnames:
                    lsel = '(%s) and (segid %s) and (resname %s)' % (sel,segid,resname)
                    mol = renumber_resid_by_resid(lsel,mol)
            else:
                lsel = '(%s) and (segid %s)' % (sel,segid)
                mol = renumber_resid_by_resid(lsel,mol)
    else:                      
        resnames = set(mol.get('resname',sel=sel))
        for resname in resnames:
            lsel = '(%s) and (resname %s)' % (sel,resname)
            mol = renumber_resid_by_resid(lsel,mol)
    return mol

def remove_apoform_ion(pdb_set,basepath = './'):
    """
    Remove ill-placed 2x50 in apoforms. 
    """

    for pdb_code in pdb_set:
        try:
            molname = glob(basepath+'receptor2curate_output/'+pdb_code+'/*apo*.pdb')[0]
            mymol = Molecule(molname, validateElements=False)
            mymol.set('name', 'O', 'resname NA')
            mymol.set('element', 'O', 'resname NA')
            mymol.set('segid', 'KO', 'resname NA')
            mymol.set('resname', 'HOH', 'resname NA')
            mymol.write(molname)
        except Exception as e:
            print("sod atom in "+pdb_code+" could not be removed because ",e)        

def covalent_ligands(mol, name, ligandsdict):
    """
    Check if the protein is covalently bound to anything similar to a ligand (not a protein, lipid or water)
    If so, check if it is a Retinol-Lysing binding. In this case change both residue names for LYR
    """
    # Take ligands marked as covalently-bound to protein residue
    covligs = [ lig for lig in ligandsdict[name] if ligandsdict[name][lig][1] ]
        
    # For every covalently-bound atom in the protein
    for lig in covligs:    
            
        # For every covalently-bound atom in the protein
        ligsel = '(resname %s)'%lig
        withinsel = '(within 2 of ('+ligsel+')) and (protein and not element H)'
        protres_resname = str(mol.get('resname', withinsel)[0])
        protres_resid = str(mol.get('resid', withinsel)[0])
        protres_chain = str(mol.get('chain', withinsel)[0])
        protres_segid = str(mol.get('segid', withinsel)[0])
        
        ligres_resid = str(mol.get('resid', 'resname '+lig)[0])
        ligres_atoms = mol.get('index', 'resname %s and resid %s'%(lig, ligres_resid))
        min_protres = min(list(mol.get('index', 'resid %s and resname %s'%(protres_resid, protres_resname))))
            
        # Reaorder atoms and resnames so ligand and residue ones are contiguous
        neworder = []
        for atom in mol.get('index', 'all'):
            if atom == min_protres:
                for my_ligres in ligres_atoms:
                    neworder.append(my_ligres)
            if atom in ligres_atoms:
                continue
            neworder.append(atom)
        mol.reorderAtoms(neworder)

        # If the ligand is a retinol and the residue a Lysine, then change both resnames to LYR
        # We have a special parameter file for this thing
        if ('RET' == lig) and (protres_resname in {'LYS','LYN'}):

            # Set resname, resid and record for new LYR residue
            mol.set('resname', 'LYR', '(resid '+protres_resid+' and resname LYS LYN) or resname RET')            
            mol.set('resid', protres_resid, 'resname LYR')
            mol.set('record', 'ATOM', 'resname LYR')
            mol.set('chain', protres_chain, 'resname LYR')
            mol.set('segid', protres_segid, 'resname LYR')

            # Change atom names to adapt to LYR nomenclature
            rettolyr = {
                'C1' : 'C17',
                'C2' : 'C16',
                'C3' : 'C15',
                'C4' : 'C14',
                'C5' : 'C12',
                'C6' : 'C11',
                'C7' : 'C10',
                'C8' : 'C9',
                'C9' : 'C80',
                'C10': 'C7',
                'C11': 'C6',
                'C12': 'C5',
                'C13': 'C3',
                'C14': 'C2',
                'C15': 'C1',
                'C16': 'C19',
                'C17': 'C18',
                'C18': 'C13',
                'C19': 'C8',
                'C20': 'C4'
            }
            # FIrst put a name+bis (to avoid re-substiutions)
            for name in mol.get('name', 'resname LYR'):
                if name in rettolyr:
                    mol.set('name', rettolyr[name]+'bis', 'name '+name)
            # ONce all substitutions have been made, replace bis
            for name in mol.get('name', 'resname LYR'):
                if 'bis' in name:
                    mol.set('name', name.replace('bis',''), 'name '+name)

        # IF not a retinol
        else:

            # Set new chain, segment, resname (COV), resid, remove Hs from protein 
            mol = renumber_covlig(mol,lig)
            mol.set('resname', 'COV', '(resid %s and resname %s) or resname %s'%(protres_resid, protres_resname, lig))
            mol.set('resid', protres_resid, 'resname COV')
            mol.set('record', 'ATOM', 'resname COV')
            mol.set('chain', protres_chain, 'resname COV')
            mol.set('segid', protres_segid, 'resname COV')

    return (mol,covligs)

def fix_and_prepare_input(inputmol,sysname,pdbcode,modresdict,isgpcr=True,first='NTER',last='CTER',prot_chain='R'):
    """
    ISMAEL FUNCTION
    Establish homogeneus nomenclature for protein residue names and segments for the system.
    """

    mol = inputmol.copy()
    aa= ' LYR LYN ALA ARG AR0 ASN ASP ASH CYX CYS CYM GLU GLH GLN GLY HIS HID HIE HIP ILE LEU LYS MET PHE PRO SER THR TRP TYR TYM VAL HSE HSD HSP '
    mol.set('resname','TIP3',sel='water')
    mol.set('chain','X',sel='resname TIP3')
    mol.set('segid','WAT',sel='water')    
    mol.set('name','OH2',sel='resname TIP3 and name OW')
    mol.set('name','H1',sel='resname TIP3 and name HW1')
    mol.set('name','H2',sel='resname TIP3 and name HW2')
    # Remove terminal Os from peptide (they give problems on capping)
    mol.remove('(resname '+aa+') and name O1 O2')
    if first == 'NTER':
        mol.set('name','HT1',sel='(resname '+aa+') and name H1')
        mol.set('name','HT2',sel='(resname '+aa+') and name H2')
        mol.set('name','HT3',sel='(resname '+aa+') and name H3')
    else:
        mol.remove('(resname '+aa+')and name H1 H2 H3')
    if last in {'CTER','CNEU','CTP','CT1'}:
        mol.set('name','OT1',sel='(resname '+aa+') and name O1')
        mol.set('name','OT2',sel='(resname '+aa+') and name O2')
        #fix
        mol.remove('(resname '+aa+') and name OT2')
    else:
        mol.set('name','O',sel='(resname '+aa+') and name O1')
        mol.remove('(resname '+aa+') and name O2')
    mol.set('name','HG1',sel='resname CYS and name HG')
    mol.set('name','HN',sel='resname HIS and name H')

    # Remove H from Retinols and Cholesterols (they will be added later during hte building process)
    mol.remove("element H and resname CLR RET")

    # Assign 'Cl' to element for Cl atoms (some chlorine atoms have 'C' as element for some reason)
    mol.set('element', 'CL', 'name CL CL1')

    his_he_resids = mol.get('resid',sel='resname HIS and name HE2')
    his_he_chains = mol.get('chain',sel='resname HIS and name HE2')
    his_he_ids = set([':'.join((chain,str(resid))) for resid,chain in zip(his_he_resids,his_he_chains)])
    his_hd_resids = mol.get('resid',sel='resname HIS and name HD1')
    his_hd_chains = mol.get('chain',sel='resname HIS and name HD1')
    his_hd_ids = set([':'.join((chain,str(resid))) for resid,chain in zip(his_hd_resids,his_hd_chains)])
    hsd_ids = his_hd_ids.difference(his_he_ids)
    hse_ids = his_he_ids.difference(his_hd_ids)
    hsp_ids = his_hd_ids.intersection(his_he_ids)
    hsd_dict = dict()
    hse_dict = dict()
    hsp_dict = dict()
    for chain,resid in [ id1.split(':') for id1 in hsd_ids]:
        if chain not in hsd_dict:
            hsd_dict[chain] = []
        hsd_dict[chain].append(resid)
    for chain,resid in [ id1.split(':') for id1 in hse_ids]:
        if chain not in hse_dict:
            hse_dict[chain] = []
        hse_dict[chain].append(resid)
    for chain,resid in [ id1.split(':') for id1 in hsp_ids]:
        if chain not in hsp_dict:
            hsp_dict[chain] = []
        hsp_dict[chain].append(resid)
    for chain in hsd_dict:
        sel1 = 'resname HIS and resid '+' '.join(hsd_dict[chain])
        mol.set('resname','HSD',sel=sel1)
    for chain in hse_dict:
        sel1 = 'resname HIS and resid '+' '.join(hse_dict[chain])
        mol.set('resname','HSE',sel=sel1)
    for chain in hsp_dict:
        sel1 = 'resname HIS and resid '+' '.join(hsp_dict[chain])
        mol.set('resname','HSP',sel=sel1)

    # Replace possible CYM residues (we dont want ionized cysteins)
    mol.set('resname', 'CYS', 'resname CYM')

    # Remove RET hydrogens if any found
    mol.remove('resname RET and element H')

    # Find modres pdbcodes in this system
    if sysname in modresdict:
        modres_string = ' '.join(modresdict[sysname])
    else:
        modres_string = ''
    # Establish segids for receptor protein
    prot_sel = '(resname '+aa+modres_string+')'
    peplig_sel = '(segid PEP L L0 L1 L2)'
    mol.set('segid', 'PX', sel=prot_sel+' and not '+peplig_sel)
    mol = autoSegment(mol,sel='segid PX', basename = 'P')
    # Establish segids for ligands (LIG for small mol and LX for peptidic)
    mol.set('segid','LIG',sel='not (water or ions or '+prot_sel+')')
    mol.set('chain','L',sel='segid LIG PEP')
    if len(mol.get('segid', 'segid PEP')):
        mol = autoSegment(mol,sel='segid PEP', basename = 'L')
    # Segids for ions and waters
    mol.set('segid','IO',sel='ions')
    mol.set('chain','N',sel='ions')

    mol = renumber_resid(mol,'water',by=1)
    mol = renumber_resid(mol,'ions',by=1)
    mol = renumber_resid(mol,'segid LIG',by=2)

    # Set 'prot_chain' as the chain id for the GPCR
    if isgpcr:
        recsegs = get_receptorsegs_frompdb(pdbcode, mol)
        mol.set('chain','T','chain '+prot_chain) # Ensure nothing else has a R chain
        mol.set('chain',prot_chain,'segid '+recsegs)

    # Delete waters within 1A of Ligand
    mol.remove('same residue as (water and within 1 of (segid LIG))')

    # Rename atoms and resnmaes of remaining waters
    isone=True
    atomsH = mol.get('index','water and element H')
    atomsO = mol.get('index', 'water and element O')
    for atom_id in atomsH:
        newname = 'H1' if isone else 'H2'
        mol.set('name', newname, 'index '+str(atom_id))
        isone = not isone
    for atom_id in atomsO:
        mol.set('name', 'OH2', 'index '+str(atom_id))
    mol.set('resname', 'TIP3', 'resname HOH')  

    protsegids = set(mol.get('segid',sel='protein'))
    return (mol,protsegids)

def make_apo(mol,recChain):
    """
    Remove all ligands or complementary proteins from system
    """
    mol.remove('not (protein or water or ion) or element Na')
    mol.remove('protein and not (chain '+recChain+')')
    prot_segids = set(mol.get('segid','protein'))
    return(mol, prot_segids)

def segchain_json(sys_mol_fixed, sysname, mystrucpath, receptor_segids_sys):
    """
    Write a file to remember to which of the original chains each segment belongs
    """
    segchain = dict()
    for seg in receptor_segids_sys:
        chain = set(sys_mol_fixed.get('chain', sel='segid '+seg))
        segchain[seg] = ' '.join(list(chain))

    with open(mystrucpath+'/segchain.json', 'w') as f:
        json.dump(segchain, f, indent = 4)
    
def chimera_superimpose(ref_pdbfile,match_pdbfile,outfile):
    """
    Superimpose a protein onto another using chimera. Only CA are taken into account
    """
    
    # Function to run Chimera commands
    def run_chimera(commands):
        chimera_path = 'chimera'  # Replace with your Chimera executable path
        cmd = f'"{chimera_path}" --nogui --script "chimera_script.py"'

        with open('chimera_script.py', 'w') as f:
            f.write('\n'.join(commands))

        os.system(cmd)

    # Define Chimera commands for superimposition
    chimera_commands = [
        'import os',
        'from chimera import runCommand as rc', 
        'from chimera import replyobj,openModels,Molecule', 
        'os.chdir(".")',
        'rc("open %s")'%ref_pdbfile,
        'rc("open %s")'%match_pdbfile,
        'rc("matchmaker #0@CA #1@CA")',
        'rc("write 1 %s")'%outfile,
        'rc("close all")'
    ]
    run_chimera(chimera_commands)
    os.remove('chimera_script.py')

def force_protonation_state(resid, target_resnames, to_resname, prep_table, mol_aligned):
    """
    Force protonation of residue with resid "resid" and resname "prev_resname" to 
    protonation state 'to_resname'. But do it only if its resname is within the 'target_resnames' set
    """
    prev_resname = mol_aligned.get('resname', sel='resid '+resid)[0]
    if prev_resname in target_resnames:
        d = prep_table
        d.loc[(d.resid == int(resid)) & (d.resname == prev_resname), 'forced_protonation'] = to_resname
    return prep_table

def find_gennum(pdbcode):
    """
    Use one of our three methods to find generic numbering
    """

    try:
        # Get standard GPCR nomenclature
        (gennum_dict,resid_dict) = find_gennum_gpcrGprot(pdbcode)
    except Exception as e:
        print("No generic numbering obtained for %s because no unrefined structure avaliable. Using alternative method... (not completely relaible!)"%(pdbcode))
        (gennum_dict,resid_dict) = find_gennum_unrefined(pdbcode)
    return(gennum_dict,resid_dict)


def find_gennum_gpcrGprot(pdbcode):
    """
    Use GPCRdb to find ResID of selected generic numbering positions in GPCRdb
    """
    
    # Download GPCRdb structure's website
    structure_data = requests.get('https://gpcrdb.org/structure/refined/'+pdbcode).content
    soup = BeautifulSoup(structure_data, 'html.parser')
    gennum_dict = {'GPCR': {}, 'Gprot':{}}
    resid_dict = {'GPCR': {}, 'Gprot':{}}

    # Extract GPCR nomenclature info
    table = soup.find('table', attrs={'id':'rotamers'})
    table_body = table.find('tbody')
    rows = table_body.find_all('tr')

    # ANalyze GPCR table with generic numbering
    for row in rows:
        cols = row.find_all('td')
        cols = [ele.text.strip() for ele in cols]
        gennum = re.sub('.\d+x','x',cols[3]) # We preffer second nomenclature type
        resid = cols[1]
        gennum_dict['GPCR'][resid] = gennum
        resid_dict['GPCR'][gennum] = resid

    # Extract GPCR nomenclature info
    table = soup.find('table', attrs={'id':'signprot_rotamers'})
    table_body = table.find('tbody')
    rows = table_body.find_all('tr')

    # ANalyze GPCR table with generic numbering
    for row in rows:
        cols = row.find_all('td')
        cols = [ele.text.strip() for ele in cols]
        gennum = cols[3]
        # gennum = re.sub('.\d+x','x',cols[3]) # We preffer second nomenclature type
        resid = cols[1]
        gennum_dict['Gprot'][resid] = gennum
        resid_dict['Gprot'][gennum] = resid

    return(gennum_dict,resid_dict)

def find_gennum_unrefined(pdbcode):
    """
    Find the Ballesteros Wanstein nomenclature for structures not yet refined in GPCRdb
    """

    # Empty dicts
    gennum_dict = {'GPCR': {}, 'Gprot':{}}
    resid_dict = {'GPCR': {}, 'Gprot':{}}
    
    # Generic numbering from this receptor via GPCRdb
    sysinfo = requests.get('https://gpcrdb.org/services/structure/'+pdbcode).json()

    # Extract name of GPCR and (if present) alpha gprot
    protname = sysinfo["protein"]
    gprot_name = False
    if ("signalling_protein" in sysinfo) and sysinfo['signalling_protein']['type']=='G protein':
        for (ent,ent_data) in sysinfo['signalling_protein']['data'].items():
            if ent_data['entry_name'].startswith('gna'):
                gprot_name =  ent_data['entry_name']

    # Resid as key, standard nomenclature thingy as value (GPCR)
    generic_nums = requests.get('https://gpcrdb.org/services/residues/extended/'+protname+'/').json()
    for pos in generic_nums:
        gennum = pos['display_generic_number']
        resid = str(pos['sequence_number'])
        if gennum:
            gennum = re.sub('.\d+x','x',gennum) # We preffer second nomenclature type
            gennum_dict['GPCR'][resid] = gennum
            resid_dict['GPCR'][gennum] = resid

    # Resid as key, standard nomenclature thingy as value (Gprot)
    if gprot_name:
        generic_nums = requests.get('https://gpcrdb.org/services/residues/extended/'+gprot_name+'/').json()
        for pos in generic_nums:
            gennum = pos['display_generic_number']
            resid = str(pos['sequence_number'])
            if gennum:
                gennum_dict['Gprot'][resid] = gennum
                resid_dict['Gprot'][gennum] = resid

    return(gennum_dict,resid_dict)

def resids_helix(pdbcode):
    """
    Filter resids belonging to helices (only one character before dot)
    """


    resids = ' '
    try:
        # Get standard GPCR nomenclature
        (gennum_dict,resid_dict) = find_gennum(pdbcode)
        sep_gen = 'x'

    except Exception as e:
        print("No refined structure found. Trying secondary method to get gennum...")
        (gennum_dict,resid_dict) = find_gennum_unrefined(pdbcode)
        sep_gen = '.'

    for (resid,gennum) in gennum_dict.items():
        gennum_split = gennum.split(sep_gen)
        if (len(gennum_split[0]) == 1) and not (gennum_split[0]=="8"):
            resids += " "+resid  

    return(resids)

def set_2x50(mol_aligned, pdbcode, thickness = None, gpcr_chain = False, sod2x50 = False):
    """
    Stablish protnation state of residue 2x50
    """
    # Download GPCRdb structure's website, and extract residue table from it
    structure_data = requests.get('https://gpcrdb.org/structure/refined/'+pdbcode).content
    soup = BeautifulSoup(structure_data, 'html.parser')
    table = soup.find('table', attrs={'id':'rotamers'})
    table_body = table.find('tbody')
    rows = table_body.find_all('tr')

    # Iterate trougth residue table to find 2.50 residue ID in this PDB file
    # Find also ids for 7.43HIS and 6.52HIS
    asp2x50 = None
    for row in rows:
        cols = row.find_all('td')
        cols = [ele.text.strip() for ele in cols]
        if cols[3] == "2.50x50":
            asp2x50 = cols[1]

    # If required, force protonation of ASP2x50 using reprepare()
    force_list = []
    asp2x50_sel = "protein and resid %s and chain %s"%(asp2x50, gpcr_chain)
    not2x50_sel = "not (%s)" % asp2x50_sel
    standard2x50 = mol_aligned.get('resname', sel=asp2x50_sel)[0]
    if sod2x50:
        force_list.append((asp2x50_sel, standard2x50))
        # prep_table = force_protonation_state(asp2x50, {'ASP','GLU'}, standard2x50, prep_table, mol_aligned)
    else:
        force_list.append((asp2x50_sel, standard2x50[:-1]+'H')) # Use protonated forms ASH and GLH 
        # prep_table = force_protonation_state(asp2x50, {'ASP','GLU'}, standard2x50[:-1]+'H', prep_table, mol_aligned)

    print(force_list)
    # Change ONLY residue 2x50, without altering anything else
    prepared_mol, prep_table = systemPrepare(mol_aligned,
                                hydrophobic_thickness=thickness,
                                no_opt = not2x50_sel,
                                no_titr = not2x50_sel,
                                no_prot = not2x50_sel,
                                return_details = True,
                                ignore_ns_errors=True,
                                force_protonation=force_list 
                                )

    # Resetting water residues to TIP3 (they are changed to TIP by proteinPrepare)
    prepared_mol.set('resname', 'TIP3', 'resname TIP')

    return(prepared_mol)

def prepare_system(mol_aligned, pdbcode, thickness = None, gpcr_chain = False, sod2x50 = False, aminergic = False, adenosine = False):
    """
    Assign protonation states using "proteinPrepare" function from HTMD, and 
    force protonation of ASP2x50 if required
    """

    if gpcr_chain:

        # Get generic numbering
        (gennum_dict,resid_dict) = find_gennum(pdbcode)

        # Iterate trougth residue table to find 2.50 residue ID in this PDB file
        # Find also ids for 7.43HIS and 6.52HIS
        gennum_gpcr = resid_dict['GPCR']
        his7x43 = gennum_gpcr['7x43'] if '7x43' in gennum_gpcr else None
        his6x52 = gennum_gpcr['6x52'] if '6x52' in gennum_gpcr else None
        asp2x50 = gennum_gpcr['2x50'] if '2x50' in gennum_gpcr else None
        asp3x32 = gennum_gpcr['3x32'] if '3x32' in gennum_gpcr else None

        # If required, force protonation of ASP2x50 using reprepare()
        force_list = []
        asp2x50_sel = "protein and resid %s and chain %s"%(asp2x50, gpcr_chain)
        standard2x50 = mol_aligned.get('resname', sel=asp2x50_sel)[0]
        if sod2x50:
            force_list.append((asp2x50_sel, standard2x50))
            # prep_table = force_protonation_state(asp2x50, {'ASP','GLU'}, standard2x50, prep_table, mol_aligned)
        else:
            force_list.append((asp2x50_sel, standard2x50[:-1]+'H')) # Use protonated forms ASH and GLH 
            # prep_table = force_protonation_state(asp2x50, {'ASP','GLU'}, standard2x50[:-1]+'H', prep_table, mol_aligned)

        # If this is an aminergic receptor, then force deprotonation 3.32 position
        asp3x32_sel = "protein and resid %s and chain %s"%(asp3x32, gpcr_chain)
        if aminergic:
            standard3x32 = mol_aligned.get('resname', sel='resid '+asp3x32)[0]
            force_list.append((asp3x32_sel, standard3x32))
            # prep_table = force_protonation_state(asp3x32, {'ASP','GLU'}, standard3x32, prep_table, mol_aligned)

        # If this is an adenosine receptor, force protonation states of HIS 6.52 and 7.43 as Hugo said
        if adenosine:
            # hiset = {'HID', 'HIS', 'HIE', 'HIP', 'HSD', 'HSP', 'HSE'}
            force_list.append(("protein and resid %s and chain %s"%(his7x43,gpcr_chain), 'HID'))
            force_list.append(("protein and resid %s and chain %s"%(his6x52,gpcr_chain), 'HIE'))
            # prep_table = force_protonation_state(his7x43, hiset, 'HID', prep_table, mol_aligned)
            # prep_table = force_protonation_state(his6x52, hiset, 'HIE', prep_table, mol_aligned)

        # Prepare protein (assign protonation states to residues)
        prepared_mol, prep_table = systemPrepare(mol_aligned,
                                    hydrophobic_thickness=thickness,
                                    return_details = True,
                                    ignore_ns_errors=True,
                                    force_protonation=force_list 
                                    )

    else:
        # Prepare protein (assign protonation states to residues)
        prepared_mol, prep_table = systemPrepare(mol_aligned,
                                    hydrophobic_thickness=thickness,
                                    return_details = True,
                                    ignore_ns_errors=True,
                                    )

    # Resetting water residues to TIP3 (they are changed to TIP by proteinPrepare)
    prepared_mol.set('resname', 'TIP3', 'resname TIP')
        
    return prepared_mol
    
def add_membrane(pdbmol,membranemol,protsegids,membrane_distance,coldist=1.3):
    # Corrections for rotational difusion
    prot = pdbmol.copy()
    protsel = 'segid '+' '.join(protsegids)
    prot.filter(protsel,_logger=False)
    r = minimalRotation(prot)
    M = rotationMatrix([0, 0, 1], r)
    pdbmol.rotateBy(M)
    pcoor = pdbmol.get('coords',sel=protsel)
    Mcoor = np.max(pcoor,axis=0)
    mcoor = np.min(pcoor,axis=0)
    # get the diagonal of XY of the protein if XY is a square 
    # which side is as long as the largest side (between X and Y) from the protein box  
    p_dim = [Mcoor[0] - mcoor[0],Mcoor[1] - mcoor[1]]
    maxXY = np.sqrt(p_dim[0]**2+p_dim[1]**2)
    minimum_box_size_x = maxXY+2
    minimum_box_size_y = minimum_box_size_x
    
    # get min max coor of the system
    minc = np.min(pdbmol.coords, axis=0).flatten()
    maxc = np.max(pdbmol.coords, axis=0).flatten()
    
    system_size_x = maxc[0] - minc[0]
    system_size_y = maxc[1] - minc[1]
    
    center_x = system_size_x/2 + minc[0]
    center_y = system_size_y/2 + minc[1]
    
    system_size = np.max([system_size_x,system_size_y])
    corr_system_size_x = np.max([minimum_box_size_x,system_size]) 
    corr_system_size_y = np.max([minimum_box_size_y,system_size])
    
    addmembdist = membrane_distance/2.0+np.max([coldist,1.5])+0.0
    
    xlim = [center_x-corr_system_size_x/2-addmembdist,center_x+corr_system_size_x/2+addmembdist]
    ylim = [center_y-corr_system_size_y/2-addmembdist,center_y+corr_system_size_y/2+addmembdist]
    
    memb = membranemol.copy()
    memb.remove('ions',_logger=False)
    memb2 = tileMembrane(memb, xlim[0], ylim[0], xlim[1], ylim[1])
    
    #from tileMembrane
    size = np.max(membranemol.get('coords', 'water'), axis=0) - np.min(membranemol.get('coords', 'water'), axis=0)
    xreps = int(np.ceil((xlim[1] - xlim[0]) / size[0]))
    yreps = int(np.ceil((ylim[1] - ylim[0]) / size[1]))
    
    membtmp_segids = ordered_unique(memb2.get('segid'))
    k=0
    for segid in membtmp_segids:
    #    memb2.set('segid','M'+str(k),sel='segid '+segid+' and not waters')
         memb2.set('segid','W'+str(k),sel='segid '+segid+' and waters')
         k += 1
            
    memb2 = renumber_segments(memb2,set(memb2.get('segid',sel='waters')),'MW')
    memb2 = renumber_segments(memb2,set(memb2.get('segid',sel='not waters')),'MEM')
    membrane_resnames = set(memb2.get('resname'))
    membrane_segids = set(memb2.get('segid'))
    
    mcenter = np.mean(memb2.get('coords',sel='segid MEM'),axis=0)
    memb2.moveBy(-mcenter)

    memb2, num = removeLipidsInProtein(pdbmol, memb2,lipidsel='lipids or waters')
    
    mol = pdbmol.copy()
    mol.append(memb2, collisions=True,coldist=coldist)
    
    return (mol, membrane_resnames,membrane_segids,xreps,yreps)

def solvate_pdbmol(mol,membrane_segids,water_thickness,water_margin,buffer=2.4,coldist=1.3,prefix='W'):
    waterbox = mol.get('coords','(waters or ions) and segid '+' '.join(membrane_segids))
    mwaterbox = np.min(waterbox, axis=0)
    Mwaterbox = np.max(waterbox, axis=0)
    coo = mol.get('coords','not (waters or ions)')
    mcoo = np.min(coo, axis=0)
    Mcoo = np.max(coo, axis=0)
    cooall = mol.get('coords','all')
    mcooall = np.min(coo, axis=0)
    Mcooall = np.max(coo, axis=0)
    #top layer
    M1z = Mcoo[2] + water_thickness/2. + np.max((coldist,buffer)) - buffer
    m1z = Mwaterbox[2] - water_margin
    M1 = [Mwaterbox[0],Mwaterbox[1],M1z]
    m1 = [mwaterbox[0],mwaterbox[1],m1z]
    print("wataerbox Max and min: ", Mwaterbox, mwaterbox)
    
    #bottom layer
    M2z = mwaterbox[2] + water_margin
    m2z = mcoo[2] - water_thickness/2.- np.max((coldist,buffer)) + buffer
    M2 = [Mwaterbox[0],Mwaterbox[1],M2z]
    m2 = [mwaterbox[0],mwaterbox[1],m2z]

    smol = solvate(mol, minmax=np.vstack((m2,M1)),prefix=prefix,buffer=buffer)

    smol, num_remove = removeAtomsInHull(smol, smol, 'name CA', 'segid "'+prefix+'[0-9]+"')
    #wtsegids = set(smol.get('segid',sel='segid "%s.*"'% prefix))
    #for segid in wtsegids:
        #smol.remove('segid %s and same resid as ( z > %g and z < %g)' % (segid,M2[2],m1[2]),_logger=False)

    return smol

def get_caps(prot_segids, mol_solvated):
    """
    Get caps for charmm build from the segment names of the protein elements 
    in the system.
    Use ACE/CT3 caps for receptor segments and NTER/CTER for non-receptor proteins (=peptide 
    and protein ligands)
    """
    caps_receptor = ['first ACE', 'last CT3']
    caps_not_receptor = ['first NTER', 'last CTER']
    nocaps = ['first none', 'last none']
    caps = dict()
    aa= ' LYR ALA ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL ASH CYM CYX GLH HID HIE HIP HSD HSE HSP LYN TYM AR0'
    for segid in prot_segids:
        caps[segid] = []
        resids = mol_solvated.get('resid','segid '+segid)
        min_res = str(min(resids))
        max_res = str(max(resids))
        # If first residue is a regular residue, cap it
        resname_min = mol_solvated.get('resname', 'segid %s and resid %s'%(segid,min_res))[0] 
        resname_max = mol_solvated.get('resname', 'segid %s and resid %s'%(segid,max_res))[0]
        if resname_min in aa:
            if segid.startswith('P'):
                caps[segid].append(caps_receptor[0])
            elif segid.startswith('L'):
                caps[segid].append(caps_not_receptor[0])
        else: 
            caps[segid].append(nocaps[1])
        # If last residue is a regular residue, cap it
        if resname_max in aa:
            if segid.startswith('P'):
                caps[segid].append(caps_receptor[1])
            elif segid.startswith('L'):
                caps[segid].append(caps_not_receptor[1])
        else: 
            caps[segid].append(nocaps[1])

    # Explicitly specify that vmd DOES NOT HAVE to cap non-protein ligands...
    ligsegids = mol_solvated.get('segid','not (protein or lipid or water or ion)')
    for segid in ligsegids:
        caps[segid] = nocaps

    return caps

def extra_parameters(pdbcode, ligandsdict, modresdict, blacklist, covligs, basepath, has_halo = False):
    """
    Get toppar files for ligand  molecules and modified residues 
    """
    
    # Make list of Ligand stringfiles (Parameters and topology)
    ligstreams = []
    for ligcode in ligandsdict[pdbcode]:
        
        #Skip blacklisted molecules and retinols or cholesterols 
        #(retinols and cholesterols already have their own parameters)
        if (ligcode not in blacklist) and not (ligcode in noparams_ligs):
            ligstreams.append('%stoppar/Ligands/%s/%s_%s_toppar.str'%(basepath,ligcode,pdbcode,ligcode))                        

    # Make list of modified-residue stringfiles (if there is any)
    modresstreams = []
    if pdbcode in modresdict:
        for modres in modresdict[pdbcode]:
            ligstreams.append('%stoppar/mod_residues/%s/%s_%s_toppar.str'%(basepath,modres,pdbcode,modres))                        

    return modresstreams+ligstreams

def add_dummy_atom(inputmol,property_dict):
    try:
        dummymol = Molecule("dummy.pdb", validateElements=False)
    except Exception: 
        dummymol = Molecule("dummy.pdb", validateElements = False)

    for prop in property_dict:
        dummymol.set(prop,property_dict[prop])
    inputmol.append(dummymol,coldist=None)
    return inputmol
def add_center_dummy_atom(inputmol,coords,property_dict):
    center = np.mean(coords,axis=0)
    property_dict['coords'] = center
    mol = add_dummy_atom(inputmol,property_dict)
    return mol
def remove_aromatic_insertions(inputmol,protsegids,coldist=1.5,outpdb=None):
    mol = inputmol.copy()
    atoms_data = [['TRP','CG CD1 CE1 NE1 CE2 CD2',5,'1'],
                 ['TRP','CD2 CE2 CZ2 CH2 CZ3 CE3',6,'2'],
                 ['PRO','N CA CB CG CD',5,''],
                 ['HIS HSD HSE HSP HID HIE HIP',
                  'CG CD1 CE1 CZ CE2 CD2 ND1 NE2',5,''],
                 ['PHE TYR TYM',
                  'CG CD1 CE1 CZ CE2 CD2',6,'']]
    beta_backup = np.copy(mol.beta)
    mol.set('beta',sequenceID((mol.resid, mol.insertion, mol.segid)))    
    
    for atom_data in atoms_data:
        atom_step = atom_data[2]
        suffix = atom_data[3]
        sel = 'resname %s and name %s' % (atom_data[0],atom_data[1])
        idxs = mol.get('index',sel=sel)
        resnames = mol.resname[idxs]
        resids = mol.resid[idxs]
        segids = mol.segid[idxs]
        coords = mol.coords[idxs,:,0]
        atom_num = len(idxs)
        if atom_num % atom_step != 0:
            raise ValueError('Missing atoms.')
        for i in range(0,atom_num,atom_step):
            property_dict = {'resname':resnames[i]+suffix,'resid':resids[i],'segid':segids[i]}
            mol = add_center_dummy_atom(mol,coords[i:i+atom_step],property_dict)
            
    if outpdb:
        mol.write(outpdb)
    
    dummy_atoms_idxs = mol.get('index',sel=dummy_sel)
    dummy_atoms_resid = mol.resid[dummy_atoms_idxs]
    dummy_atoms_segid = mol.segid[dummy_atoms_idxs]
    removed_indexes = []
    for idx,resid,segid in zip(dummy_atoms_idxs,dummy_atoms_resid,dummy_atoms_segid):
        r_idx1 = mol.get('index', sel='not (%s) and same beta as ((exwithin %g of (index %d)) and not (resid %d and segid %s) and not protein)'  % (dummy_sel,coldist,idx,resid,segid))
        removed_indexes = removed_indexes + r_idx1.tolist()
        
    if len(removed_indexes) > 0:
        removed_indexes_str = ' '.join(str(x) for x in removed_indexes)
        mol.remove('index '+removed_indexes_str)
    mol.remove(dummy_sel,_logger=False)
    inv_idx1 = np.setdiff1d(np.arange(len(beta_backup)), np.array(removed_indexes), assume_unique=True)
    mol.beta = beta_backup[inv_idx1]
        
    print('WARNING: removed '+str(len(removed_indexes))+' atoms within '+str(coldist)+' of a protein aromatic ring')
    
    return (mol,removed_indexes)

def equilwrap_structure(equildir):
    """
    Build PDB structure from last frame of equillibrated structure
    """
    # Skip if already done
    outfile = equildir+'equillibrated.pdb'
    if os.path.exists(outfile):
        return
    
    mol = Molecule(equildir+'structure.psf', validateElements=False)
    mol.read(equildir+'structure.pdb')
    mol.read(equildir+'output.xtc')
    gpcr_sel = "protein and chain P"

    # Remove all frames except last
    copmol = mol.copy()
    copmol.dropFrames(keep=int(mol.numFrames-1))
    
    # Wrap system
    copmol.wrap(gpcr_sel)

    # Align frames
    copmol.align('all', refmol=Molecule(equildir+'structure.pdb'))
    
    copmol.write(outfile)

################# This function is old and hasnt been used in a while
################# Consider updating if you ever need it back
# def mutate(mol, pdbcode, equildir, mutdir, mutations, protchain, basepath, topparpath):
#     """
#     Create mutant version of structure with PDBcode
#     """

#     # Skip if mutant already exists
#     if os.path.exists(mutdir+'structure.pdb'):
#         return
    
#     # Load and mutate structure according to what the dict says
#     selection_mutated = "chain %s and resid "%protchain
#     for mut in mutations:
#         mol.mutateResidue('chain %s and resid %s'%(protchain,mut[0]), mut[1])
#         selection_mutated += mut[0]+' '

#     # Use proteinPrepare to put side-chains in mutated residues
#     # proteinPrepare will be applied only to protein, ligand and protein waters
#     nonprot_chains_set = set(mol.get('chain', sel='not protein'))
#     nonprot_chains = 'chain '+' '.join(item for item in nonprot_chains) 
#     scafmol = mol.copy()
#     mol.remove(nonprot_chains)
#     scafmol.filter(nonprot_chains)
#     # prepared_mol = systemPrepare(mol, ## Temporarely disabled because acellera people cant do anything stright
#     #                       no_opt= "not (%s)" %selection_mutated,
#     #                       no_titr= "not (%s)" %selection_mutated,
#     #                       no_prot= "not (%s)" %selection_mutated,
#     #                      )
#     prepared_mol = mol
#     prepared_mol.append(scafmol)
#     prepared_mol.set('resname', 'TIP3', 'resname TIP')

#     #Remove 5 CLA and SOD atoms to leave charmm.build some margin to redo system charge
#     cla_remove = ' '.join(map(str, prepared_mol.get('resid',"segid I and resname CLA")[:5]))
#     sod_remove = ' '.join(map(str, prepared_mol.get('resid',"segid I and resname SOD")[:5]))        
#     prepared_mol.remove("segid I and resid "+cla_remove+sod_remove)

#     # Check if system has lone-pair hallogen atoms. If it does, use legacy CGenFF parameters
#     (cgenff_par,cgenff_top,has_halo) =get_params_cgenff(mol, topparpath)

#     #Obtain extra parameters for ligands and modified residues 
#     ligstreams=extra_parameters(pdbcode, ligandsdict, modresdict, blacklist, [], basepath, has_halo)

#     # Caps
#     caps = get_caps(set(mol.get('segid', 'protein and chain %s'%protchain)), mol)

#     #Build system with changes
#     buildmol = charmm.build(prepared_mol, 
#                                topo=topos+cgenff_top, 
#                                param=params+cgenff_par,
#                                stream=streams+ligstreams,
#                                outdir=mutdir,
#                                caps=caps,
#                                saltconc=0.15)

def define_equilibration(const_sel = 'protein and name C CA N O or not (protein or lipid or water or ions) and noh or segid IO WAT and noh', 
                         simtime = 40, minimize = 5000, timestep = 2, temperature = 310):
    """
    Default restraint selection includes main-chain atoms, crystal waters/ions and ligand molecules
    """

    restrs = [AtomRestraint(const_sel, 2, [(1,"0"),(1,"%dns" % int(simtime*0.5)),(0,"%dns" % int(simtime*0.75))], "xyz")]
    md = Equilibration()
    md.restraints = restrs
    md.runtime = simtime
    md.timeunits = 'ns'
    md.temperature = temperature
    # md.nvtsteps = 0
    md.acemd.barostatconstratio = 'on'
    md.acemd.minimize = minimize
    md.acemd.restart = 'off'
    md.acemd.timestep = timestep
    return md

def define_production(timestep, trajperiod, temperature=310, runtime=500):
    md = Production()
    md.runtime = runtime
    md.timeunits = 'ns'
    md.temperature = temperature
    md.acemd.timestep = timestep
    md.acemd.barostatconstratio = 'on'
    md.acemd.restart = 'off'
    md.acemd.trajectoryperiod = trajperiod
    md.acemd.bincoordinates = 'output.coor'
    md.acemd.extendedsystem  = 'output.xsc'
    md.acemd.binvelocities = 'output.vel'
    return md

def job_commands(sourcedir, nodedir):
    """
    Commands to copy files into dwarf, run ACEMD there, move everything back to /gpcr/
    And then exit the bash script
    """
    return [
            'cd',
            'mkdir -p '+nodedir,
            'mv %s* %s'%(sourcedir, nodedir),
            'touch %ssimrunning' %sourcedir,
            'eval "$(/opt/miniconda3/bin/conda shell.bash hook)"',
            'eval "$(conda shell.bash hook)"',
            'conda activate acemd4_mod',
            'cd '+nodedir,
            'bash run.sh',
            'cd ..',
            'mv %s* %s'%(nodedir, sourcedir),
            'rm -r %s %ssimrunning'% (nodedir,sourcedir),
            'exit'
            ]

def transmem_atoms(strucpath):
    """
    Extract Z-coordinate of P atoms of phospholipids, make the averages of top and bottom layers, 
    select the CA inbetween and return VMD selection of them 
    """
    mol = Molecule(strucpath, validateElements=False)
    z_coors = mol.get('coords',sel='lipid and name P')[:, 2]
    z_coor_av = np.average(z_coors)
    top_z_av = np.average([ z for z in z_coors if z>z_coor_av ])
    bot_z_av = np.average([ z for z in z_coors if z<z_coor_av ])
    transmem_serial = mol.get('serial',sel='name CA and protein and (z >= %s) and (z <= %s)'%(bot_z_av, top_z_av))
    transmem_serial_str = ' '.join([str(a) for a in transmem_serial]) 
    transmem_sel = 'serial '+transmem_serial_str
    return(transmem_sel)


def wrap_alig(strucpdb, mymol_pdb, proddir, prot_sel, transmem_sel):
    """
    Wrap and align the system around the transmembrane region
    """

    # To avoid repeating wrapping in Trajectories already wrapped, check the existance of this file
    outname = proddir+'output_wrapped.xtc'
    if os.path.exists(outname):
        print('replicate already wrapped. Skipping...')
        return
    
    # Wrap system
    mymol = Molecule('%sstructure.psf' % (proddir), validateElements=False)
    mymol.read(proddir+'output.xtc')
    mymol.wrap(prot_sel)

    # Align frames and write
    mymol.align(transmem_sel, refmol=mymol_pdb)
    mymol.write(outname)


def wrap_alig_vmd(topopath, strucpath, trajpath, outname, prot_sel, transmem_sel):
    """
    Very big systems cannot be properly wrapped with MDanalysis (it core-dumps).
    FOr those cases, we need to do some shenanigands with VMD
    """
    outdcd = outname+'.dcd'
    outxtc = outname+'.xtc'
    cmd = ['vmd', topopath]+[' -dispdev', 'text']
    vmd = Popen(cmd, stdin=PIPE, universal_newlines=True)
    vmd.stdin.write("\n".join([
    ' package require pbctools',
    'set molid top',
    'mol addfile %s type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all'%trajpath,
    '# Function to wrap',
    'pbc wrap -center com -centersel "%s" -compound residue -all'%transmem_sel,
    'pbc wrap -center com -centersel "%s" -compound residue -all'%transmem_sel,
    'pbc wrap -center com -centersel "%s" -compound residue -all'%prot_sel,
    'pbc wrap -center com -centersel "%s" -compound residue -all'%prot_sel,
    '# Function to align',
    'mol addfile %s type pdb first 0 last 1 step 1 filebonds 1 autobonds 1 waitfor all'%strucpath,
    'set transmem_sel "%s"'%transmem_sel,
    'set ref [atomselect 0 $transmem_sel frame last]',
    'set sel [atomselect 0 $transmem_sel]',
    'set all [atomselect 0 all]',
    'set n [molinfo top get numframes]',
    'for { set i 0 } { $i < $n } { incr i } {',
    '$sel frame $i',
    '$all frame $i',    
    '$all move [measure fit $sel $ref]',
    '}',
    'set orig_lastframe [ expr $n-2 ]',
    '# Save trajectory in dcd (silly vmd cannot save xtc) and a PDB of the first frame (for GPCRmd)',
    'animate write dcd {%s} beg 0 end $orig_lastframe skip 0 waitfor all'%outdcd,
    'quit',
    ]))
    vmd.stdin.close()
    vmd.wait()

    # Convert dcd to xtc and remove original dcd
    os.system("mdconvert %s -f -o %s; rm %s"%(outdcd,outxtc,outdcd))


################# Single-use functions
def remove_spare_hydrogens(pdb_set,basepath = './'):
    """
    Remove protonation Hydrogens from unprotonated residue names in curator's structures
    """
    for pdb_code in pdb_set:
        for pdbfile in glob(basepath+'receptor2curate_output/'+pdb_code+'/*.pdb'):
            try:
                mymol = Molecule(pdbfile, validateElements=False)
                mymol.remove( 'resname ASP and name HD2')
                mymol.remove( 'resname GLU and name HE2')
                mymol.remove( 'resname HID and name HE2')
                mymol.remove( 'resname HIE and name HD1')
                mymol.remove( 'resname CYM and name HG1')
                mymol.remove( 'resname LYN and name HZ3')
                mymol.remove( 'resname TYM and name HH')
                mymol.remove( 'resname AR0 paramand name HH12')
                mymol.write(pdbfile)

            except Exception as e:
                print("sod atom in "+pdb_code+" could not be removed because ",e)    

def count_atoms(pdb_set, basepath):
    """
    Prints the total number of atoms in each build system
    """
    for a in pdb_set:
        try:
            mol = Molecule(basepath+'simulation_output/build/'+a+'_apo/structure.pdb', validateElements=False)
            l = len(mol.get('resid','all'))
            print(a,' '+str(l))

        except Exception as e:
            pass    

def moecurated_paramchem(username,password,ligpath,ligcode,pdbcode):
    """
    Submit ligand to paramchem to obtain topology-parametrs file
    Version for ligands manually protonated with MOE
    """

    # Make two paramchems: one with legacy parameters and the other with modern ones (for hallogen systems)
    for c_option in ['c','']:
        
        #Pattern to find non-HTMD-compatible lines
        lph_pat = re.compile('^ATOM.*LPH|LONEPAIR')

        # Ligand mol2file
        mol2file = ligpath+pdbcode+"_protonlig.mol2"
        # topology-parametetrs file output
        topparfile_path = ligpath+pdbcode+"_legacy_toppar.str" if c_option else ligpath+pdbcode+"_toppar.str"
        # Errors and warnings file output
        erwar_path = ligpath+"legacy_paramchem_stder.txt" if c_option else ligpath+"legacy_paramchem_stder.txt"
    
        # Omit ligand if its toppar file already exists
        if os.path.exists(topparfile_path):
            leg = 'legacy' if c_option else 'latest'
            print('%s toppar for ligand %s already exists. Skipping...' % (leg,ligcode))
            continue
        
        # Ligands for which chimera creates wrong topology have to be autodetermined in paramchem
        badtop_lig = {'7AB'}
        b_option = 'b' if ligcode in badtop_lig else ''

        #Open ligandfile to upload in paramchem
        myfile = {
                'filename': open(mol2file)
        }

        # Define POST variables for login in Paramchem
        datalogin = {
            'usrName': username,
            'curPwd': password,
            'rndNo': str(random.randrange(100000000, 999999999)),
            'submitBtn': 'Submit',
        }

        # Variables for paramchem
        myparam = {
                #'param_a': 'a', #Include parameters usually included in newer versions of CGenff (versions that we cant use)
                'param_b' : b_option , # Make paramchem determine bond order (active if covalent-bound ligand)
                'param_c': c_option, # Use CGenFF legacy v1.0, for HTMD is not yet prepared for newer Charmm versions             
        }

        # Open web session
        with requests.Session() as s:

            # Login into paramchem
            response_login = s.post('https://cgenff.silcsbio.com/userAccount/userWelcome.php', 
                       data=datalogin,
                       verify=False)

            # Submit our ligand molecule into paramchem
            response_upload = s.post('https://cgenff.silcsbio.com/initguess/processdata.php', 
                       files = myfile,
                       data = myparam,
                        )

            # Download Topology-parameters of our molecule file from paramchem, and store it.
            # But remember submissions in paramchem are limited weekly
            # Download also stderr, just in case
            match = re.search('<path>(\w+)</path>', response_upload.text)
            if match:
                code = match.group(1)
                response_ligfile = s.get('https://cgenff.silcsbio.com/initguess/filedownload.php?file='+code+'/'+pdbcode+'_protonlig.str')
                response_stder = s.get('https://cgenff.silcsbio.com /initguess/filedownload.php?file='+code+'/'+pdbcode+'_protonlig.err')
                with open(topparfile_path, 'wb') as topparfile:
                        topparfile.write(response_ligfile.content)
                with open(erwar_path, 'wb') as erwar:
                        erwar.write(response_stder.content)                            
            else:
                print(response_upload, response_upload.content)
                print('Your paramchem account has reached its weekly submission limit. Please, intrudce a new account or wait to the next monday to continue')            
                return

            #Delete lines with LPH (new feature from CHARMM not tolerated by HTMD)
            """
            with open(topparfile_path, "r") as f:
                lines = f.readlines()
            with open(topparfile_path, "w") as f:
                for line in lines:
                    if not re.match(lph_pat, line):
                        f.write(line)
            """
###################################
## Part 5: Upload results to GPCRmd
###################################

#### Functions
def resp_to_dict(resp):
    # Convert a json reponse into a dictionary
    return eval(resp.content.decode('UTF-8').replace('true', 'True').replace('false','False'))

def get_headers(s, subm_id):
    """
    Obtain headers, csrftoken and sessionid from current session object
    """
    sessionid = str(s.cookies['sessionid'])
    csrftoken = str(s.cookies['csrftoken'])
    headers = {
        'Referer' : mainurl+'/dynadb/step1/'+subm_id,
        'Cookie' : 'csrftoken=%s; sessionid=%s' %(csrftoken,sessionid),
        'X-CSRFToken' : csrftoken,
    }
    return(sessionid, csrftoken, headers)

def find_dyn_name(protdict,pdbcode,ligdict,apo):
    """
    Find an appropiate name for this submitted system
    """

    # Find GPCR name (if any)
    gpcr_name = ''
    for (chain,protdata) in protdict.items():
        if protdata['prot_type'] == 1:
            gpcr_name = protdata['uniname']

    # Get first ligand molecule at random 
    lig_name = ''
    if len(ligdict.keys()):
        first_lig_resname = next(iter(ligdict.items()))[0]
        first_name = ligdict[first_lig_resname]['name']
        first_inchikey = ligdict[first_lig_resname]['inchikey']
        lig_name = first_name if len(first_name)<29 else first_inchikey 
    
    # Stablish system name
    if gpcr_name:
        if apo:
            dynname = gpcr_name+' (apoform)'
        elif lig_name:
            dynname = gpcr_name+' in complex with '+lig_name
        else:
            dynname = gpcr_name
    else:
        apoword = 'apoform' if apo else 'complex'
        dynname = pdbcode+apoword
    
    # It doesnt accept more than 100 characters
    dynname = dynname[0:97]+'...'

    return(dynname)

def check_chains(pdbcode, mymol):
    """
    Check how many chains from the original PDB structure remain in mymol structure
    And to which Segments of our structure they correspond
    """
    # Load blosum score matrix to align proteins
    blosum62 = substitution_matrices.load("BLOSUM62")

    # Obtain sequences for original PDB file chains, and classifying them by chainID
    # Also check which segment corresponds to what chain
    pdbmol = load_mol_pdbcode(pdbcode)
    chainseg = {}
    chainset = set(pdbmol.get('chain', sel='protein'))
    pdbmol_segseqs = pdbmol.sequence()
    pdbmol_chainseqs = {}
    for chain in chainset:
        segid = np.unique(pdbmol.get('segid', sel='protein and chain '+chain))[0]
        pdbmol_chainseqs[chain] = pdbmol_segseqs[segid]

    # Merging all protein chains in our systems into a single megachain
    mymol_megachain = ''
    for seg,chain in mymol.sequence().items():
        mymol_megachain = mymol_megachain + chain

    # Aligning sequences from original PDB to simulated PDB megasequence
    # To know which chains from the original PDB are preserved
    chain_present = {}
    segtochain = {}
    # For each chain in the original PDB file
    for chain,seq in pdbmol_chainseqs.items():
        chain_present[chain] = False
        pdb_len = len(seq)
        # For each segment in our molecule
        for seg,myseq in mymol.sequence().items():

            # Skip alignment if there is a huge size difference between both sequences
            mylen = len(myseq)
            lenratio = mylen/pdb_len 
            if (lenratio>5) or lenratio<0.2:
                continue
            
            # Align our molecule segments to the pdb chains
            aligs = pairwise2.align.localms(seq, myseq, 5,-1, -1.5, 0)
            # Check if any of the alignments has more than 4 score by position
            shortest_len = mylen if mylen<pdb_len else pdb_len
            is_present = any([ (alig[2]/shortest_len) > 4 for alig in aligs ]) 
            if is_present:
                chain_present[chain] = is_present
                segtochain[seg] = chain
                break

    return (chain_present, segtochain)

def get_uniprot_alignment(mymol, segpos, uniprotkbac, isoform, entry_segs):
    """
    Find out the exact coordinates of each Uniprot segment in our protein
    """

    # Obtain sequence of uniprotkbac
    uniprot_id = "%s-%s"%(uniprotkbac,isoform) if isoform else uniprotkbac 
    response = requests.get("https://rest.uniprot.org/uniprotkb/"+uniprotkbac+".fasta")
    uniseq = ''.join(response.text.split('\n')[1:-1])
    # Match this uniprot to our protein segments, and save matches
    for seg,myseq in mymol.sequence().items():

        # Align our molecule segments to the pdb chains
        aligs = pairwise2.align.localms(myseq, uniseq, 5,-1, -1, -0.5)
        # For every alignment
        for alig in aligs:
            start = alig.start
            end = alig.end
            # If the alignment has no significant amount of mismatches, keep it as segment
            if (alig.score/(end-start)) > 3.5:
                
                # Find coordinates of this segment
                alig_myseq = alig[0]
                
                # Find start of uniprot sequence in alginment
                if alig_myseq.startswith('-'):
                    starting = 0
                else:
                    starting = start - 1 
                # Find ending of uniprot sequence in alginment            
                if alig_myseq.endswith('-'):            
                    ending = len(myseq) - 1
                else:
                    ending = end - 1 - alig.seqA.count('-')
                
                if uniprotkbac not in entry_segs:
                    entry_segs[uniprotkbac] = []
                    
                # Do not add repeated alignments                
                beg = segpos[seg][starting]
                end = segpos[seg][ending]
                repeated_seg = any([ ((seg['beg']==beg) and (seg['end']==end)) for seg in entry_segs[uniprot_id]])
                if not len(entry_segs[uniprot_id]) or not repeated_seg:
                    entry_segs[uniprot_id].append({
                                'chainid' : mymol.get('chain', 'segid '+seg)[0],
                                'segid' : seg,
                                'beg' : segpos[seg][starting],
                                'end' : segpos[seg][ending]})

    return(entry_segs) 

def get_segids_chain(mol, chainid):
    """
    Extract segids of a certain chain in a molecule htmd object
    """
    segids = set(mol.get('segid',sel='chain '+chainid))
    result = ' '.join(str(item) for item in segids)
    return(result)

def get_receptorsegs_frompdb(pdbcode, mymol):
    """
    Get informarion of the system specified in the pdbcode from the RCSB-PDB webpage.
    In this one only protein information, nothing about ligands
    """

    # Get residue names in our system
    allresnames = set(mymol.resname)
    
    # Extract information from pdb webpage using api
    ligdict = dict()
    protdict = dict()
    datadict = dict()
    entry_segs = dict()
    segpos = dict()

    # Obtain dictionary of sorted protein resids in each segment
    for seg in set(mymol.get('segid', 'protein')):
        segpos[seg] = list(dict.fromkeys(mymol.get('resid', 'segid '+seg)))
    
    # Get information from PDB api
    pdbdict = requests.get('https://data.rcsb.org/graphql?query={\
        entry(entry_id: "'+pdbcode+'") {\
            polymer_entities {\
                entity_poly{pdbx_strand_id}\
                rcsb_polymer_entity{pdbx_description}\
                rcsb_polymer_entity_align{\
                    aligned_regions {length}\
                    reference_database_accession \
                    reference_database_isoform \
                }\
            }\
        }\
    }').json()['data']
                    
    # Extract protein chains information
    receptor_segs=""
    for poly in pdbdict['entry']['polymer_entities']:

        # Determine if this polymer (chain) is a GPCR
        uniname = poly['rcsb_polymer_entity']['pdbx_description'].lower()
        chainIds = poly['entity_poly']['pdbx_strand_id'].split(',')
        isgpcr = any([(name in uniname) for name in gpcr_names])
        # Check list of uniprots that are GPCR but have no "gpcr-name" (e.g.: uncharacterized protein)
        if poly['rcsb_polymer_entity_align']:
            for alreg in poly['rcsb_polymer_entity_align']:
                uniprotkbac = alreg['reference_database_accession']
                if any([(name == uniprotkbac) for name in gpcr_uniprots]):
                    isgpcr = True
        
        # Find segments of Receptor protein
        if isgpcr:
            # Get uniprots present in the system, and their corresponbding segments
            for alreg in poly['rcsb_polymer_entity_align']:
                uniprotkbac = alreg['reference_database_accession']
                isoform = alreg['reference_database_isoform']
                uniprot_id = "%s-%s"%(uniprotkbac,isoform) if isoform else uniprotkbac 
                entry_segs = get_uniprot_alignment(mymol, segpos, uniprotkbac, isoform, entry_segs)
                # If there was an alignment match
                if uniprotkbac in entry_segs:
                    for segchain in entry_segs[uniprotkbac]:
                        receptor_segs += segchain['segid']+' '

    return(receptor_segs)

def get_pdb_info(pdbcode, mymol, ligandsdict):
    """
    Get informarion of the system specified in the pdbcode from the RCSB-PDB webpage.
    Mainly uniprot sequences, chainIDs and uniprot codes for the chains
    """
    
    # Make string for selecting lignads in rcsb's api
    ligs = '"'+'","'.join(ligandsdict[pdbcode].keys())+'"'
    
    #Check which chains from the original PDB structure are preserved in mymol structure
    (chain_present, segtochain) = check_chains(pdbcode, mymol)
    
    # Get ligand molecules present in our system (basically anything that is not a protein)
    ligset = set(mymol.get('resname', sel='not protein'))
    
    # Extract information from pdb webpage using api
    ligdict = dict()
    protdict = dict()
    datadict = dict()
    
    # Get information from PDB api
    pdbdict = requests.get('https://data.rcsb.org/graphql?query={\
        entry(entry_id: "'+pdbcode+'") {\
            exptl {method}\
            polymer_entities {\
                entity_poly{pdbx_strand_id}\
                rcsb_polymer_entity{pdbx_description}\
                rcsb_polymer_entity_align{\
                    aligned_regions {length}\
                    reference_database_accession \
                    reference_database_isoform \
                }\
            }\
        }\
        chem_comps(comp_ids: ['+ligs+']) {\
            rcsb_chem_comp_descriptor{InChIKey}\
            rcsb_chem_comp_descriptor{InChI}\
            rcsb_chem_comp_descriptor{SMILES} \
            chem_comp{id,name,pdbx_formal_charge} \
        }\
    }').json()['data']
    
    # Get ligand information and classify it (in case the system actually has ligands)
    if len(pdbdict['chem_comps'])>0:
        for lig in pdbdict['chem_comps']:
            InchiKey = lig['rcsb_chem_comp_descriptor']['InChIKey']
            Inchi = lig['rcsb_chem_comp_descriptor']['InChI']
            smiles = lig['rcsb_chem_comp_descriptor']['SMILES']
            Resname = lig['chem_comp']['id']
            Name = lig['chem_comp']['name']
            Charge = lig['chem_comp']['pdbx_formal_charge']
    
            # Exclude ligands not present in our simulated system
            if Resname in ligset:
                ligdict[Resname] = {
                    'inchikey' : InchiKey,
                    'inchi' : Inchi,
                    'smiles' : smiles,
                    'name' : Name,
                    'charge' : Charge
                }

    # Extract protein chains information
    for poly in pdbdict['entry']['polymer_entities']:

        # Skip if no alignment with uniprot entity avaliable
        if not poly['rcsb_polymer_entity_align']:
            continue
    
        # Get alignment information in order to determine which uniprot segment we should take
        prevallen = 0
        uniprotkbac = ''
        isoform = ''
        for alreg in poly['rcsb_polymer_entity_align']:
            allen =  sum([seg['length'] for seg in alreg['aligned_regions']])
            unicode = alreg['reference_database_accession']
            myisoform = alreg['reference_database_isoform']
            if allen > prevallen:
                uniprotkbac = unicode
                isoform = myisoform
            else:
                uniprotkbac = uniprotkbac
                isoform = isoform
            prevallen = allen            
            
        uniname = poly['rcsb_polymer_entity']['pdbx_description'].lower()
        chainIds = poly['entity_poly']['pdbx_strand_id'].split(',')

        # Determine protein type of this chain
        if (('ceptor' in uniname) or ('opsin' in uniname) or ('frizzled' in uniname) or ('smoothen' in uniname)):
            prot_type = 1 # GPCR
        elif ('Guanine nucleotide-binding protein' in uniname):
            prot_type = 2 # G-protein subunit
        elif ('Beta-arrestin' in uniname):
            prot_type = 3 # Beta arrestin
        else:
            prot_type = 4 # We'lll assume any other protein as a peptide ligand 
        
        # If this chain is present in our mymol structure
        for chainId in chainIds:
            if chain_present[chainId]:
                # Get which segment(s) this chain is assigned to. SKip if none
                segs = [ seg for seg,chain in segtochain.items() if chain == chainId ]
                if not segs:
                    continue
                # Order data for this chain
                protdict[chainId] = {
                    'uniname' : uniname,
                    'uniprotkbac' : uniprotkbac,
                    'isoform' : isoform if isoform else '1',
                    'prot_type' : prot_type,
                    'segs' : segs,
                }

    # Determine experimental method used, and use the corresponding id in GPCRmd database
    method = pdbdict['entry']['exptl'][0]['method'].lower()
    if 'x-ray' in method: 
        method_id = 0
    elif 'nmr' in method:
        method_id = 1
    elif 'electron microscopy' == method:
        method_id = 4
    else:
        method_id = 5 # Other method
            
    return (protdict,ligdict,segtochain,method_id)

def login(s, username='XXXXXXXXXXXX', password="XXXXXXXXXXXX"):

    if (username=='XXXXXXXXXXXX') or (password=="XXXXXXXXXXXX"):
        raise Exception("please define arguments 'username' and 'password' for loging into GPCRmd")

    # Login into GPCRmd 
    # Will you go, lassie, go? And we'll all gooo toghetheeeer
    headers = {
        'Cookie': 'csrftoken=JzdqGFC5HJ2Gk4gxksv5IjFv2FNeUXKW',
        'Referer': mainurl+'/accounts/login/',
    }
    datalogin = {
        'username': username,
        'password': password,
        'next' : '/accounts/memberpage/',
        'csrfmiddlewaretoken' : 'JzdqGFC5HJ2Gk4gxksv5IjFv2FNeUXKW'
    }
    print('loging into GPCRmd')
    logo = s.post(mainurl+'/accounts/login/', 
               data=datalogin,
               headers=headers)
    
    return(s)

def new_submission(s, mainurl):
    """
    Create a new submission entry in GPCRmd
    """
    response_new = s.get(mainurl + '/dynadb/step0/')
    soup = BeautifulSoup(response_new.text, 'html.parser')
    step1_link = soup.find('form',attrs={'id' : 'mainform'}).get('action')
    subm_id = step1_link.split('/')[-2]
    print('new submission %s created'%subm_id)
    return subm_id

def check_completeness(prodpath):
    """
    Check if we actually have the three necessary replicates for submission
    """

    incomplete = False
    for trajid in range(1,repnum+1):
        trajpath = '%s/rep_%d/output_wrapped.xtc'%(prodpath,trajid)
        if not os.path.exists(trajpath):
            incomplete = True

    return(incomplete)

def new_step1(s, subm_id, dynname, pdbcode, trajperiod, timestep, prodpath, method_id, apo):
    """
    Fullfill and send the step 1 of the new submissionform
    # Looong liiie the fieeeelds of Athenryyyyy...
    # Where once, we watch, the small seabirds flyyyyy....
    """
    # Get headers and stuff
    (sessionid, csrftoken, headers)=get_headers(s, subm_id)
    
    # All the stuff we need to send 
    apoword = 'apoform' if apo else 'complex'
    sysname = pdbcode+' '+apoword
    data_submit = {
        'name' : dynname,
        'type' : 0 if apo else 1, # 0 for apo, 1 for complex
        'pdbid' : pdbcode,
        'description' : 'Classical unbiased (NVT ensemble) complex flexibility assay.',
        'source_type' : method_id,
        'state_type' : 0, # We'll implement this later
        'id_dynamics_methods': '1', #Molecular mechanics,
        'software': 'ACEMD3',
        'sversion': 'GPUGRID', # TODO: Is that correct??,
        'ff': 'CHARMM',
        'ffversion': 'c36 Jul 2021',
        'id_assay_types': 1, # Orthosteric (un)/binding,
        'id_dynamics_membrane_types': 2, # Homogeneus membrane,
        'id_dynamics_solvent_types': 2, # TIP3P solvent,
        'timestep': timestep,
        'delta': (trajperiod*timestep)/1e6,
        'add_info': 'autosubmission', # time/frames
        'gpcrcom' : True, # Is this gpcr part of the gpcrmd community effort (yes, the whole 3th round is)
        'passcode' : 'weareGPCRs',
    }
    files_submit = {
        'dynamics' : open(prodpath+'structure.pdb', 'rb'),
    }

    # Submit step 1
    rep = s.post(mainurl+'/dynadb/step1_submit/'+subm_id+'/',
            headers = headers,
            files = files_submit,
            data = data_submit)
    print('step1 finalized ',rep)

    
def smalmol_getdata(s, molid, inchikey, subm_id, moltype, resname, sdfpath="", pdbdata=""):
    """
    Get data for small molecules 
    """
    
    # Get headers and stuff
    (sessionid, csrftoken, headers)=get_headers(s, subm_id)

    # Get info of molecule
    rep = s.get(mainurl+'/dynadb/smalmol_info/?inchikey=%s&submission_id=%s'%(inchikey,subm_id),
        headers = headers,
        )
    rep_dict = rep.json()
    

    # Prepare submission dictionary with the data of this molecule
    if rep_dict['smiles']:
        mol_datasubmit = {
            "smalmol_type"+molid : moltype,
            "smalmol_inchikey"+molid : inchikey,
            "smalmol_name"+molid : rep_dict['name'],
            "smalmol_resname"+molid : resname,
            "smalmol_iupac"+molid : rep_dict['iupac'],
            "smalmol_chemblid"+molid :  rep_dict['chemblid'],
            "smalmol_cid"+molid :  rep_dict['cid'],
            "smalmol_inchi"+molid : rep_dict['inchi'],
            "smalmol_sinchi"+molid :  rep_dict['sinchi'],
            "smalmol_sinchikey"+molid :  rep_dict['sinchikey'],
            "smalmol_smiles"+molid :  rep_dict['smiles'],
            "smalmol_netcharge"+molid :  rep_dict['net_charge'],
            "smalmol_synonims"+molid :  rep_dict['other_names'],
            "smalmol_description"+molid :  '',
            "image_path"+molid :  rep_dict['imagepath']
        }
    # If no response, use data from pdb
    else: 
        mol_datasubmit = {
            "smalmol_type"+molid : moltype,
            "smalmol_inchikey"+molid : inchikey,
            "smalmol_name"+molid : pdbdata['name'][:60],
            "smalmol_resname"+molid : resname,
            "smalmol_iupac"+molid : pdbdata['name'],
            "smalmol_chemblid"+molid :  '',
            "smalmol_cid"+molid :  '',
            "smalmol_inchi"+molid :pdbdata['inchi'],
            "smalmol_sinchi"+molid : pdbdata['inchi'],
            "smalmol_sinchikey"+molid : pdbdata['inchikey'],
            "smalmol_smiles"+molid : pdbdata['smiles'],
            "smalmol_netcharge"+molid : pdbdata['charge'],
            "smalmol_synonims"+molid :  '',
            "smalmol_description"+molid :  '',
            "image_path"+molid : ''
        }

    # Add SDFfile to the submission if required
    mol_files = {} if rep_dict['inGPCRmd'] else {"sdfmol"+molid : open(sdfpath)} 
    
    return(mol_datasubmit,mol_files)
    
def new_step2(s, subm_id, ligdict, basepath='./'):
    """
    Fullfil and send step2 of new submission form
    """    
    print('initiating step 2: small molecule data')

    # Get headers and stuff
    (sessionid, csrftoken, headers)=get_headers(s, subm_id)
    
    # Introduce common small molecules (waters, lipids and ions)
    i = 0
    num_entries = []
    step2_data = {
        'csrfmiddlewaretoken' : csrftoken,
        'submit' : 'submit'
    }
    step2_files = {}
    common_mols = {
        'TIP3' : ('XLYOFNOQVPJJNP-UHFFFAOYSA-N',6),
        'POPC' : ('WTJKGGKOPKCXLL-VYOBOKEXSA-N',7),
        'SOD' : ('FKNQFGJONOIPTF-UHFFFAOYSA-N',8),
        'CLA' : ('VEXZGXHMUGYJMC-UHFFFAOYSA-M',8)
    }
    sdfpath = basepath+'smalmol_sdfs/'
    for resname,(inchikey,moltype) in common_mols.items():
        sdfile = sdfpath+resname+'.sdf'
        (mol_datasubmit, mol_files) = smalmol_getdata(s, str(i), inchikey, subm_id, moltype, resname, sdfile)
        step2_data.update(mol_datasubmit)
        step2_files.update(mol_files)
        num_entries.append(i)
        i+=1

    # Introduce ligands: all will be defined as orthosteric ligands
    lignames = []
    for lig,ligdata in ligdict.items():

        # Extract basic info
        name = ligdata['name']
        inchikey = ligdata['inchikey']

        # Avoid blacklisted molecules or ions
        if (lig in detergent_blacklist) or (lig in glucids_blacklist)  or (len(lig) == 2):
            continue
        # Download ligand, store it into temporary file
        response = requests.get('https://files.rcsb.org/ligands/view/'+lig+'_ideal.sdf')
        with open('tmpfile.sdf','wb') as tmpout:
            tmpout.write(response.content)

        # CLR (cholesterol) is usually a experimental lipid 
        moltype = '3' if lig=='CLR' else '0'
            
        # Retrieve molecule information
        (mol_datasubmit, mol_files) = smalmol_getdata(s, str(i), inchikey, subm_id, moltype, lig, 'tmpfile.sdf', ligdata)
        step2_data.update(mol_datasubmit)
        step2_files.update(mol_files)
        num_entries.append(i)
        i+=1
        os.remove('tmpfile.sdf')
    
    # Add number of entries
    step2_data['num_entries'] = ','.join(map(str, num_entries))
    
    # Replace empty entries with an empty string
    for key,val in step2_data.items():
        if not val:
            step2_data[key] = ''
    
    # Submit step2
    rep = s.post(mainurl+'/dynadb/step2_submit/'+subm_id+'/',
        headers = headers,
        data = step2_data,
        files = step2_files)
    print('step2 finalized ',rep)

   
def new_step3(s, subm_id, protdict):
    """
    Fullfil and sent the step3 (protein chains) of the new submission form
    """
    print('initiating step 3: protein chains')
    
    # Get headers and stuff
    (sessionid, csrftoken, headers)=get_headers(s, subm_id)
    
    # Retrieve segment information using GPCRmd submission form
    entrydict = s.get(mainurl+'/dynadb/find_prots/'+subm_id+'/',
        headers = headers).json()
    
    segments = {seg['segid'] : seg for entry in entrydict for seg in entry['segments']}
    
    # Start dict with the data to be send into GPCRmd
    # But thougth they paved the footways here with goldust, I still would chose... my Isle of Innesfree
    step3_data = {
        'csrfmiddlewaretoken':csrftoken,
        'submission_id' : subm_id,
        'submit' : 'Submit',
    }
    
    # Iterate over pdb chains. For each one we'll create an entry
    i=0
    entrynum = []
    for chain,chain_data in protdict.items():
        ii = str(i)
        uniprotkbac = chain_data['uniprotkbac']
        isoform = chain_data['isoform']
        chainsegs = chain_data['segs']
    
        # Retrieve info using uniprot code
        unidict = s.get(mainurl+'/dynadb/prot_info/'+subm_id+'/',
                          params = {
                              'uniprotkbac' : uniprotkbac,
                              'isoform' : isoform,
                              'submission_id' : subm_id,
                          },
                          headers = headers
                         ).json()
    
        # We'll obtain the alignment of our protein sequence and its UNIPROT sequence using
        # GPCRmd's tool for do it
        data_alig = {
            'unisequence' : unidict['Sequence'],
            'submission_id' : subm_id
        }
    
        # Get segments for alignment data
        j = 0
        segnums = []
        for seg in chainsegs:
            seginfo = segments[seg]
            ss = str(j)
            data_alig['chain'+ss]= seginfo['chain']
            data_alig['segid'+ss]= seg
            data_alig['from'+ss]= seginfo['resid_from']
            data_alig['to'+ss]= seginfo['resid_to']
            segnums.append(ss)
    
            # Store segment information as well in the submit dict
            step3_segdata = {
                ss+'seg'+ii: j,
                ss+'pdbid'+ii: seginfo['pdbid'],
                ss+'sourcetype'+ii: seginfo['source_type'],
                ss+'chain'+ii: seginfo['chain'],
                ss+'segid'+ii: seg,
                ss+'from_resid'+ii: seginfo['resid_from'],
                ss+'to_resid'+ii: seginfo['resid_to'],
            }
            step3_data.update(step3_segdata)
    
            j += 1
        data_alig['segnums']=','.join(segnums)
        step3_data['num_segs'+ii]=','.join(segnums)
    
        # Send alignment request
        alignment = s.post(mainurl+'/dynadb/get_alignment/',
                      data = data_alig,
                      headers = headers
                     ).json()
    
        # Get mutations
        mutations_dict = s.post(mainurl+'/dynadb/protein/get_mutations/',
                      data = {
                          'alignment' : alignment['alig_fa'],
                          'sequence' : unidict['Sequence'],
                      },
                      headers = headers
                     ).json()['mutations']
    
        # Store mutations in submission dictionary
        h = 0
        mutnums = []
        for mut in mutations_dict:
            hh = str(h)
            step3_mutdata = {
                hh+'from'+ii : mut['from'],
                hh+'to'+ii : mut['to'],
                hh+'resid'+ii : mut['resid'],
            }
            step3_data.update(step3_mutdata)
            mutnums.append(hh)
            h+=1
        step3_data['num_muts'+ii]=','.join(mutnums)
    
        # Get chain name according to uniprot
        names_array = unidict['Protein names'].split(" (")
        names_array = list(map(lambda x: re.sub("\)$", "", x), names_array)) # Convert the name string into name aray
        name = names_array.pop(0) # We get the first name as the official one
        other_names = ';'.join(names_array)
        
        # Start storing information to send in step3
        step3_entrydata = {
            'isoform'+ii : isoform,
            'prot_uniprot'+ii : uniprotkbac,
            'prot_type'+ii : unidict['prot_type'],
            'name'+ii : name,
            'aliases'+ii : other_names,
            'species_code'+ii : unidict['Organism (ID)'],
            'species_name'+ii : unidict['Organism'],
            'unisequence'+ii : unidict['Sequence'],
        }
        step3_data.update(step3_entrydata)
    
        entrynum.append(ii)    
        i+=1
    
    # Store number of entries
    step3_data['num_entries'] = ','.join(entrynum)
    
    # Submit step3
    rep = s.post(mainurl+'/dynadb/step3_submit/'+subm_id+'/',
        headers = headers,
        data = step3_data)
    
    print('step3 finalized ',rep)

def new_step4(s, subm_id, prodpath, repath):
    """
    Fullfil and sent the step4 (files) of the new submission form
    """
    print('initiating step 4: simulation files')

    # Get headers and stuff
    (sessionid, csrftoken, headers)=get_headers(s, subm_id)

    # Prepare trajectories
    step4_files = []
    for trajid in range(1,repnum+1):
        trajpath = '%s/rep_%d/output_wrapped.xtc'%(prodpath,trajid)
        step4_files.append(('trj', ('traj%d.xtc'%trajid, open(trajpath, 'rb'))))
        # trajpath = '%s/rep_%d/trajtest.dcd'%(prodpath,trajid)
        # step4_files.append(('trj', ('trajtest.dcd', open(trajpath, 'rb'))))
    
    # Prepare the rest of files (parameters files will name-changed to avoid problems with filextensions)
    step4_files.append(
        ('dyn', ('structure.psf', open(repath+'structure.psf','rb')))
    )
    step4_files.append(
        ('prm', ('parameters.prm', open(repath+'parameters','rb')))
    )
        
    # Encode data
    m=MultipartEncoder(fields=step4_files)
    headers['Content-Type'] = m.content_type
    
    # Submit xxzzzstep4
    rep = s.post(mainurl+'/dynadb/step4_submit/'+subm_id+'/',
        headers = headers,
        data = m)
    
    print('step4 finalized ',rep)

def new_step5(s, subm_id, pub_code, provo_pmid):
    """
    Fullfil and sent the step4 (files) of the new submission form
    """
    # Get headers and stuff
    (sessionid, csrftoken, headers)=get_headers(s, subm_id)

    data_submit = {
        'pub_code' : pub_code,
        'pmid' : provo_pmid,
    }
    
    # Submit step 1
    rep = s.post(mainurl+'/dynadb/step5_submit/'+subm_id+'/',
            headers = headers,
            data = data_submit)
    print('step5 finalized ',rep)

def create_subm_ids_file(pdb_set,resultspath,basepath):
    """
    Create json file with subm_ids of submitted systems
    """
    subm_ids = {}
    for pdbcode in pdb_set:
        for apo in ['','_apo']:
            sysname = pdbcode+apo
            submfile = "%s/production/%s/submitted.txt"%(resultspath,sysname)
            if os.path.exists(submfile):
                with open(submfile,'r') as infile:
                    subm_id = infile.readline()
                subm_ids[sysname] = subm_id

    with open(basepath+'gandalf_subm_ids.json','w') as outf:
        json.dump(subm_ids, outf, indent=4)

def get_dynid(subm_id):
    """
    Obtain dyn_id from subm_id of a simulation
    """
    resp = requests.get('https://gpcrmd.org/api/search_sub/2253').json()
    return(resp[0]['dyn_id'])