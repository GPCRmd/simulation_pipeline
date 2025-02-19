# Load functions, modules and global variables required in our pipeline
import requests
from bs4 import BeautifulSoup
from importlib import reload  
from sys import argv
from getpass import getpass
import os
import re
import tarfile
import argparse as ap
import json
import traceback
from requests_toolbelt import MultipartEncoder

############
## Functions
############

def json_dict(path):
	"""
	Converts json file to pyhton dict.
	"""
	json_file=open(path)
	json_str = json_file.read()
	json_data = json.loads(json_str)
	return json_data

def get_headers(s, subm_id):
    """
    Obtain headers, csrftoken and sessionid from current session object
    """
    sessionid = str(s.cookies['sessionid'])
    csrftoken = str(s.cookies['csrftoken'])
    headers = {
        'Referer' : mainurl+'/dynadb/step1/'+subm_id,
        'Cookie' : 'csrftoken=%s; sessionid=%s' %(csrftoken,sessionid),
        'X-CSRFToken' : csrftoken
    }
    return(sessionid, csrftoken, headers)

def login(s, username='XXXXXXXXXXXX', password="XXXXXXXXXXXX"):

    if (username=='XXXXXXXXXXXX') or (password=="XXXXXXXXXXXX"):
        raise Exception("please define arguments 'username' and 'password' for loging into GPCRmd")

    # Login into GPCRmd 
    # Will you go, lassie, go? And we'll all gooo toghetheeeer
    headers = {
        'Cookie': 'csrftoken=II5zdSvPpkVzLV3siBFHv5kGdtbsZr2C',
        'Referer': mainurl+'/accounts/login/',
    }
    datalogin = {
        'username': username,
        'password': password,
        'next' : '/accounts/memberpage/',
        'csrfmiddlewaretoken' : 'II5zdSvPpkVzLV3siBFHv5kGdtbsZr2C'
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

def autosubmit_step1(s, subm_id, sys_info, pdbpath):
	"""
	Fullfill and send the step 1 of the new submissionform
	# Looone liiie the fieeeelds of Athenryyyyy...
	# Where once, we watched, the small seabirds flyyyyy....
	"""
	# Get headers and stuff
	(sessionid, csrftoken, headers)=get_headers(s, subm_id)
	
	# Dropdown options
	systypes = ["APO","COMPLEX"]
	srctypes = ["XRAY","NMR","EM","OTH"]
	dynmethods = ["","MM","QM","NMR","MONTECARLO","MSM","METADYN","UMBRELLA","BIAS","SCALED"]
	assaytypes = ["","ORTO","ACT","OLIGO","ALLOBIND","ALLOMOD","GPROTCOUP","BARRCOUP"]
	memtypes = ["","IMP","HOMO","HETE"]
	solvtypes = ["","IMP","TIP3","TIP4","TIPS","SPC","SPCE","OTH"]

	# All the stuff we need to send 
	data_submit = {
		'name' : sys_info['name'],
		'type' : systypes.index(sys_info['type']),
		'pdbid' : sys_info['pdbid'],
		'description' : sys_info['description'],
		'source_type' : srctypes.index(sys_info['source_type']),
		'id_dynamics_methods': dynmethods.index(sys_info['id_dynamics_methods']),
		'software': sys_info['software'],
		'sversion': sys_info['sversion'],
		'ff': sys_info['ff'],
		'ffversion': sys_info['ffversion'],
		'id_assay_types': assaytypes.index(sys_info['id_assay_types']),
		'id_dynamics_membrane_types': memtypes.index(sys_info['id_dynamics_membrane_types']),
		'id_dynamics_solvent_types': solvtypes.index(sys_info['id_dynamics_solvent_types']),
		'timestep': sys_info['timestep'],
		'delta': sys_info['delta'],
		'add_info': sys_info['add_info'],
	}
	files_submit = {
		'dynamics' : open(pdbpath, 'rb')
	}

	# Submit step 1
	rep = s.post(mainurl+'/dynadb/step1_submit/'+subm_id+'/',
			headers = headers,
			files = files_submit,
			data = data_submit)
	print('step1 finalized ',rep)

	# Raise exceptions if something failed
	if not rep.ok:
		s.post(mainurl+'/dynadb/delete_submission/'+subm_id)
		print(traceback.format_exc())
		raise Exception('Step 1 submitted returned %d'%rep.status_code)

def smalmol_getdata(s, molid, inchikey, subm_id, moltype, resname, sdfpath=""):
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
    
    # Download ligand, store it into temporary file, and load it into submit if necessary
    if rep_dict['inGPCRmd']:
        mol_files = {}
    else:
        response = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/SDF"%rep_dict['cid'])
        with open('tmpfile.sdf','wb') as tmpout:
            tmpout.write(response.content)
        mol_files = {"sdfmol"+molid : open('tmpfile.sdf')} 
    
    return(mol_datasubmit,mol_files)

def autosubmit_step2(s, subm_id, smalmols):
	"""
	Fullfil and send step2 of new submission form
	"""    
	print('initiating step 2: small molecule data')

	# Get headers and stuff
	(sessionid, csrftoken, headers)=get_headers(s, subm_id)
	smalmoltypes = ["ORT","ALLO","EION","ELIP","EWAT","OTHC","BWAT","BLIP","BION","BOTH"]

	# Introduce common small molecules (waters, lipids and ions)
	i = 0
	num_entries = []
	step2_data = {
		'csrfmiddlewaretoken' : csrftoken,
		'submit' : 'submit'
	}
	step2_files = {}
	# Introduce ligands: all will be defined as orthosteric ligands
	lignames = []
	for smalmol in smalmols:
		   
		# Retrieve molecule information
		moltype = smalmoltypes.index(smalmol['type'])
		(mol_datasubmit, mol_files) = smalmol_getdata(s, str(i), smalmol['inchikey'], subm_id, moltype, smalmol['resname'])
		step2_data.update(mol_datasubmit)
		step2_files.update(mol_files)
		num_entries.append(i)
		i+=1
	
	# Add number of entries
	step2_data['num_entries'] = ','.join(map(str, num_entries))
	
	# Replace empty entries with an empty string
	for key,val in step2_data.items():
		if (not val) and (val!= 0):
			step2_data[key] = ''
	
	# Submit step2
	rep = s.post(mainurl+'/dynadb/step2_submit/'+subm_id+'/',
		headers = headers,
		data = step2_data,
		files = step2_files)
	print('step2 finalized ',rep)

	# Raise exceptions if something failed
	if not rep.ok:
		s.post(mainurl+'/dynadb/delete_submission/'+subm_id)
		print(traceback.format_exc())
		raise Exception('Step 2 submitted returned %d'%rep.status_code)


def autosubmit_step3(s, subm_id, pdbcode, protlist):
	"""
	Fullfil and sent the step3 (protein chains) of the new submission form
	"""
	print('initiating step 3: protein chains')

	# Get headers and stuff
	(sessionid, csrftoken, headers)=get_headers(s, subm_id)
	segsrc = ["XRAY","NMR","ABIN","HOMO","THR","MD","EM","OTH"]
	prot_types = ["GPCR","GPROT","BARR","PEP","OTH"]
	
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
	for prot in protlist:
		ii = str(i)
		uniprotkbac = prot['uniprotkbac']
		isoform = prot['isoform']
		segs = prot['segs']
		prot_type_name = prot['prot_type']

		# Retrieve info using uniprot code
		# If no uniprot code retrieved from input json
		if uniprotkbac:
			unidict = s.get(mainurl+'/dynadb/prot_info/'+subm_id+'/',
							  params = {
								  'uniprotkbac' : uniprotkbac,
								  'isoform' : isoform,
								  'submission_id' : subm_id,
							  },
							  headers = headers
							 ).json()
			# Get chain name according to uniprot
			names_array = unidict['Protein names'].split(" (")
			names_array = list(map(lambda x: re.sub("\)$|\) \[.*\]$", "", x), names_array)) # Convert the name string into name aray
			name = names_array.pop(0) # We get the first name as the official one
			other_names = ';'.join(names_array)            
		else: 
			unidict = {
			'Organism (ID)' : prot['organism_id'],
			'Organism' : prot['organism'],
			'Sequence' : '',
			'prot_type' : prot_types.index(prot_type_name),
			}
			name = prot['name']
			other_names = prot['other_names']

		# We'll obtain the alignment of our protein sequence and its UNIPROT sequence using
		# GPCRmd's tool for do it
		data_alig = {
			'unisequence' : unidict['Sequence'],
			'submission_id' : subm_id
		}

		# Get segments for alignment data
		j = 0
		segnums = []
		for seg in segs:
			ss = str(j)
			data_alig['chain'+ss]= seg['chainid']
			data_alig['segid'+ss]= seg['segid']
			data_alig['from'+ss]= seg['beg']
			data_alig['to'+ss]= seg['end']
			segnums.append(ss)

			# Store segment information as well in the submit dict
			step3_segdata = {
				ss+'seg'+ii: j,
				ss+'pdbid'+ii: pdbcode,
				ss+'sourcetype'+ii: segsrc.index(seg['source_type']),
				ss+'chain'+ii: seg['chainid'],
				ss+'segid'+ii: seg['segid'],
				ss+'from_resid'+ii: seg['beg'],
				ss+'to_resid'+ii: seg['end'],
			}
			step3_data.update(step3_segdata)

			j += 1
		data_alig['segnums']=','.join(segnums)
		step3_data['num_segs'+ii]=','.join(segnums)

		# Alignment and mutations thing only works if there is an uniprot
		if uniprotkbac:

			# Send alignment request
			alignment = s.post(mainurl+'/dynadb/get_alignment/',
						  data = data_alig,
						  headers = headers
						 ).json()['alig_fa']

			# Get mutations
			mutations_dict = s.post(mainurl+'/dynadb/protein/get_mutations/',
						  data = {
							  'alignment' : alignment,
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

		
		# Start storing information to send in step3
		step3_entrydata = {
			'isoform'+ii : isoform,
			'prot_uniprot'+ii : uniprotkbac,
			'name'+ii : name,
			'aliases'+ii : other_names,
			'species_code'+ii : unidict['Organism (ID)'],
			'species_name'+ii : unidict['Organism'],
			'unisequence'+ii : unidict['Sequence'],
			'prot_type'+ii : unidict['prot_type'],
		}
		step3_data.update(step3_entrydata)

		entrynum.append(ii)    
		i+=1

	# Store number of entries
	step3_data['num_entries'] = ','.join(entrynum)

	# Submit step3
	print(step3_data)
	rep = s.post(mainurl+'/dynadb/step3_submit/'+subm_id+'/',
		headers = headers,
		data = step3_data)
	
	print('step3 finalized ',rep)

	# Raise exceptions if something failed
	if not rep.ok:
		#s.post(mainurl+'/dynadb/delete_submission/'+subm_id)
		print(traceback.format_exc())
		raise Exception('Step 3 submitted returned %d'%rep.status_code)


def autosubmit_step4(s, subm_id, files):
	"""
	Fullfil and sent the step4 (files) of the new submission form
	"""
	print('initiating step 4: simulation files')

	# Get headers and stuff
	(sessionid, csrftoken, headers)=get_headers(s, subm_id)

	# Prepare trajectories
	step4_files = []
	for trajfile in files['trajs']:
		step4_files.append(('trj', (os.path.basename(trajfile), open(trajfile, 'rb'))))
		#trajpath = '%s/rep_%d/trajtest.dcd'%(prodpath,trajid)
		#step4_files.append(('trj', ('trajtest.dcd', open(trajpath, 'rb'))))
	
	# Prepare the rest of files (parameters files will name-changed to avoid problems with filextensions)
	step4_files.append(
		('dyn', (os.path.basename(files['topo']), open(files['topo'],'rb')))
	)
	step4_files.append(
		('prm', (os.path.basename(files['params']), open(files['params'],'rb')))
	)
	# step4_files.append(
		# ('prt', (os.path.basename(files['proto']), open(files['proto'],'rb')))
	# )
	# Other files are optional, so upload only if there is file
	if ('other' in files) and (files['other']) :
		step4_files.append(
			('oth', (os.path.basename(files['other']), open(files['other'],'rb')))
		)

		# Encode data
	m=MultipartEncoder(fields=step4_files)
	headers['Content-Type'] = m.content_type
	
	# Submit step4
	rep = s.post(mainurl+'/dynadb/step4_submit/'+subm_id+'/',
		headers = headers,
		data = m)
	
	print('step4 finalized ',rep)

	# Raise exceptions if something failed
	if not rep.ok:
		s.post(mainurl+'/dynadb/delete_submission/'+subm_id)
		print(traceback.format_exc())
		raise Exception('Step 4 submitted returned %d'%rep.status_code)


def autosubmit_step5(s, subm_id, ref_info):
	"""
	Fullfill and send the step 1 of the new submissionform
	# Ouur loove was on the wind
	# We had dreaaams and songs to sing
	# It's so lonely, down the fields, of Athenry
	"""

	# Get headers and stuff
	(sessionid, csrftoken, headers)=get_headers(s, subm_id)
	
	# All the stuff we need to send 
	data_submit = {
		'doi' : ref_info['doi'],
		'authors' : ref_info['authors'],
		'title' : ref_info['title'],
		'pmid' : ref_info['pmid'],
		'journal' : ref_info['journal'],
		'issue' : ref_info['issue'],
		'volume' : ref_info['volume'],
		'pages' : ref_info['pages'],
		'pub_year' : ref_info['pub_year'],
		'url' : ref_info['url'],

	}
	files_submit = {
		'dynamics' : open(pdbpath, 'rb')
	}

	# Submit step 1
	rep = s.post(mainurl+'/dynadb/step1_submit/'+subm_id+'/',
			headers = headers,
			files = files_submit,
			data = data_submit)
	print('step1 finalized ',rep)


##############################
## Initial variables and stuff
##############################

# Arguments
parser = ap.ArgumentParser(description="this calculates interaction frequencies for given simulation")
parser.add_argument(
	'-i',
	dest='input_json',
	action='store',
	required=True,
	default=False,
	help='Json file with information regarding your simulations to submit. You can find a template in input_subm.json'
)
parser.add_argument(
	'-u',
	dest='username',
	action='store',
	required=True,
	default=False,
	help='Username of your GPCRmd account.'
)
parser.add_argument(
	'--url',
	dest='url',
	action='store',
	default='https://submission.gpcrmd.org/',
	help='Database to submit your simulations into (in case David is testing stuff in his localhost).'
)
args = parser.parse_args()
input_json = args.input_json
mainurl = args.url
username = args.username

# Ask user to introduce GPCRmd password
print('please intrudouce your GPCRmd password')
password = getpass()

# Load data
input_dict = json_dict(input_json)

## Step -2: Login into GPCRmd
with requests.Session() as s:
	login(s, username, password)

# For each of the currently-working-with systems defined in Part 1
# for pdbcode in pdb_set:
for entry in input_dict:    

	try:

		sysname = entry['sys_info']['name']
		print('\n############### Submitting '+sysname+' simulation ########################')            
		
		## Step 0: Create a new submission
		subm_id = new_submission(s, mainurl)
		# subm_id = "1477"
		print('submission id: '+subm_id)

		## Step 1: General information
		autosubmit_step1(s, subm_id, entry['sys_info'], entry['files']['coords'])

		## Step 2: Small molecules 
		autosubmit_step2(s, subm_id, entry['smallmolecules'])

		## Step 3: Protein chains
		autosubmit_step3(s, subm_id, entry['sys_info']['pdbid'], entry['proteins'])

		## Step 4: Dynamics information
		autosubmit_step4(s, subm_id, entry['files'])
	
	except Exception as e:
		print("Simulation "+sysname+" could not be submitted because ",e)
		print(traceback.format_exc())
