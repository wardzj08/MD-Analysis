import pandas as pd
import numpy as np


mu3H = list(pd.read_csv('mu3HParticles.txt', delim_whitespace = True, header = None)[0])
mu3O = list(pd.read_csv('mu3OParticles.txt', delim_whitespace = True, header = None)[0])
print(mu3H)
print(mu3O)

#skip 73 lines on read in
# rewrite pparameters for O and H manually, including adding the new mu3H parameter
# Then use this to rewrite individual atoms
# if mu3H, rewrite as new seperate H atom type (atom 7) with new charge parameter (0.435)
# if mu3O, keep type, but change charge? (to -0.7000)
# Leave bond parameters as are for now
fr = 'data/data.UIO-66-2.ddec'
fw = 'data.UiO-66-mu3MOD.ddec'
with open(fr, 'r') as r:
	lines = r.readlines()
	#print(lines[73])
	with open(fw, 'w') as w:
		# First 74 lines remain the same (pair coefficients modified manually)
		for line in lines[:74]:
			w.write(line)
		end = 74 + 3648
		# Modify the individual atoms (change atom type and charges as needed)
		for line in lines[74:end]:
			#print(int(line.split()[0]))
			# Change mu3 Hydrogen values
			if int(line.split()[0]) in list(mu3H):
				print(line.split()[0])
				typeH = 6 # Use this to add new hydrogen atom type for only mu3H, resulting in 2 different hydrogen types
				chargeH = 0.43500 # Use to alter the charge on these hydrogens  (to match IPA trappe)
				splt = line.split()
				# Write new line
				w.write(f'     {splt[0]}      {splt[1]}        {typeH}     {chargeH}      {splt[4]}      {splt[5]}      {splt[6]}\n')
				continue
			# Alter mu3-O atoms. Change charge on atoms as needed. Atom type stays the same
			if int(line.split()[0]) in mu3O:
			#	print(line.split()[0])
				chargeO = -0.70000
				line = line.split()
				w.write(f'     {line[0]}      {line[1]}        {line[2]}     {chargeO}      {line[4]}      {line[5]}      {line[6]}\n')
				continue
			elif int(line.split()[0]) not in mu3O or int(line.split()[0]) not in mu3H:
				#print(line.split()):q

				# Any other atoms can just be rewritten as is
				w.write(line)

		# All lines after the atom list do not need to be changed (bonds, angles, dihedrals stay the same)
		for line in lines[end:]:
			w.write(line)