import os
import re

# Open file with atom codes of all the lacking parameters
infile = open('./warnings.txt', 'r')


dic_das = {
	'CT' : 'CG301', # No-hydrogens alyphatic carbon (C)
	'CT1' : 'CG311', # One-hydrogen alyphatic carbon (CH)
	'CT2' : 'CG321', # Two-hydrogens alpyphatic carbon (CH2)
	'CT2A' : 'CG321', # Same as before (GLU, HSP chi1/chi2 fitting)
	'CT3' :'CG331', # Three-hydrogens alyphatic carbon (methyl group) (CH3)
	'C' : 'CG2O1', # Amide carbon (the one in peptidic bonds)
	'NH1' : 'NG2S1', # Peptide nitrogen
	'NH3' : 'NG3P3', # NH3+ terminal amino nitrogen
	'H' : 'HGP1', # Polar H
	'O' : 'OG2D1', # Amide oxygen (the one in peptidic bonds)
	'HB1' : 'HGA1', # Aliphatic H in alkane (CH)
	'HB2' : 'HGA2', # Aliphatic H in alkane (CH2)
	'HA1' : 'HGA1', # Aliphatic H in alkane (CH) (new LJ parameters)
	'HA2' : 'HGA2', # Aliphatic H in alkane (CH2) (new LJ parameters)
	'CP1' : 'CG3C51', # Proline CA
	'CP2' : 'CG3C52', # Proline CB/CG
	'CP3' : 'CG3C52', # Proline CD
	'N' : 'NG2S0', # Proline N
	'NP' : 'NG3P2' # Proline N-terminal NH2+
}

outlines = {'BONDS' : {}, 'ANGLES' : {}, 'DIHEDRALS' : {}, 'IMPROPERS' : {}}
numpat = re.compile('^\w+\s+\w+\s+\w+\s+\w+\s+(\d+)')
warnpat = re.compile('(for|improper|dihedral)(.*)')

# For each lacking parameter named in the file
for line in infile:

	# Omit crossterm lines (have to be prepared separatedly)
	if 'crossterm' in line:
		continue

	# Remove warning and leave olny atom names
	line = re.sub('^#.*(for|improper|dihedral)','', line)
	line = line.rstrip()+' '
	line = line.upper()
	line_len = len(line.split())
	greplines = set()
	repl_words = set()
	word_counter = dict()
	nword_list = []
	for word in line.split():
		if word in dic_das:
			nword = dic_das[word]
			repl_words.add(word)
		else:
			nword = word
		nword_list.append(nword)
		word_counter[nword] = word_counter[nword]+1 if nword in word_counter else  1

	# Create grep commnad line to search the parameters of this atoms
	grepmatch = 'grep -ih "%s.* \! \|%s.* \! " \
		*/toppar.str \
		*/legacy_toppar.str \
		../Ligands/*/????_toppar.str \
		../TOP_PARAMS_ACE3/General_top_params/parameters/par_all36_cgenff.prm | sort'%('.*'.join(nword_list),'.*'.join(nword_list[::-1]))
	greps = os.popen(grepmatch+'| sort').read().split('\n')
	if len(greps) == 1:
		print(line,greps)
		pass

	# Iterate over the grep line results
	linepat = '(^'+('\w+\s+'*line_len)+')(\d+.*)!'
	lpat = re.compile(linepat)
	for grepline in greps:
		# If is the parameter type we are looking for 
		match = lpat.search(grepline)
		if match:
			str_atoms = match.group(1)
			values = match.group(2)
			str_atomsl = str_atoms.split()
			str_atoms = '  '.join(str_atomsl)
			for word in repl_words:
				ntimes = line.count(word+' ')
				str_atoms = str_atoms.replace(dic_das[word], word, ntimes)
			
			# Create reverse-atoms version of our parameter line (e.g.: a B C A -> A B C a)
			rev_atoms = match.group(1)
			rev_atomsl = rev_atoms.split()
			rev_atomsl.reverse()
			rev_atoms = '  '.join(rev_atomsl)
			for word in repl_words:
				ntimes = line.count(word+' ')
				rev_atoms = rev_atoms.replace(dic_das[word], word, ntimes)
			rev_atomsl = rev_atoms.split('  ')
			rev_atomsl.reverse()
			rev_atoms = '  '.join(rev_atomsl)			

			for atoms in [rev_atoms, str_atoms]:
				if line_len == 2:
					if atoms in outlines['BONDS']:
						outlines['BONDS'][atoms].add(values)
					else:
						outlines['BONDS'][atoms] = {values}
				elif line_len == 3:
					if atoms in outlines['ANGLES']:
						outlines['ANGLES'][atoms].add(values)
					else:
						outlines['ANGLES'][atoms] = {values}
				elif line_len == 4:
					strength = numpat.search(grepline).group(1)
					if int(strength) >10: # UNREALIBLE
						if atoms in outlines['IMPROPERS']:
							outlines['IMPROPERS'][atoms].add(values)
						else:
							outlines['IMPROPERS'][atoms] = {values}
					else:
						if atoms in outlines['DIHEDRALS']:
							outlines['DIHEDRALS'][atoms].add(values)
						else:
							outlines['DIHEDRALS'][atoms] = {values}

# Forcing some more parameters in the dictionary
outlines['IMPROPERS']['N  CG2O1  CP1  CP3'] = {'0.0000         0      0.0000 ! manually added'}
outlines['DIHEDRALS']['NH3  CT1  C  NG2S1'] = {'0.4000  1     0.00 ! manually added'}

# Write new parameters file
outfile = open('/gpcr/users/daranda/doctorat/GPCR_simulations/toppar/TOP_PARAMS_ACE3/David_top_params/parameters/modres_params.prm', 'w')
#outfile = open('modres_params.prm','w')
outfile.write("""
read param card flex append
* Hand-made parameters David Aranda from Paramchem data, to make compatible modified residue's topologies-parameters with normal residue's ones
* CHARMM General Force Field (CGenFF) program version 2.4.0
*

	""")

for key in outlines:
	outfile.write("\n"+key+"\n")
	for atoms in outlines[key]:
		for values in outlines[key][atoms]:
			outfile.write(atoms+'  '+values+"\n")

outfile.write("""

END
RETURN
""")

outfile.close()

"""
#WARNING: illegal value: no parameters exist for dihedral cg2o1 nh1 ct2 c
#WARNING: illegal value: no parameters exist for dihedral ct1 c ng2s1 cg301
#WARNING: illegal value: no parameters exist for dihedral ct1 ct1 c ng2s1
CG2O1  CG311  NG2S1  C      0.2000  1   180.00 ! PROT ala dipeptide, new C VDW Rmin, adm jr., 3/3/93c
C  CT1  NH1  CG2O1      0.2000  1   180.00 ! PROT ala dipeptide, new C VDW Rmin, adm jr., 3/3/93c
CG2O1  CG301  NG2S1  C      0.2000  1   180.00 ! 9DT , from CG2O1 CG311 NG2S1 CG2O1, PENALTY= 8
"""