import pandas as pd
import numpy as np


Remove_types = [10, 11, 12, 13]
Remove_bonds = [10, 11, 12]
Remove_angles = [13,14,15]
Remove_dihedrals = [11]


# Remove
fr = 'potentials/UiO-66-20ipa222-Q.data'
fw = 'potentials/RoggeUiO66Framework.data'
reorder = True
if reorder:
	curr_atom = 1
	num_atoms = 4348

	curr_bond = 1
	num_bonds = 5296

	curr_angle = 1
	num_angle = 7936

	curr_di = 1
	num_di = 7856
	oldAtomNum = []
	newAtomNum = []

with open(fr, 'r') as r:
	lines = r.readlines()
	#print(lines[73])
	with open(fw, 'w') as w:
		# First 74 lines remain the same (pair coefficients modified manually)
		for line in lines[:100]:
			w.write(line)
		#end = 74 + 3648
		# Modify the individual atoms (change atom type and charges as needed)
		for line in lines[100:4552]:
			# Remove adsorbate atoms
			print(line.split())
			if len(line.split()) <= 2:
				w.write(line)
				continue
			if  int(line.split()[2]) in Remove_types:
				continue

			elif int(line.split()[2]) not in Remove_types:
				# Any other atoms can just be rewritten as is
				line = line.split()
				oldAtomNum += [int(line[0])]
				newAtomNum += [curr_atom]
				w.write(f'     {curr_atom}      {line[1]}        {line[2]}     {line[3]}      {line[4]}      {line[5]}      {line[6]}\n')
				curr_atom += 1

		for line in lines[4552:9932]:
			# Remove adsorbate atoms
				print(line.split())
				if len(line.split()) <= 2:
					w.write(line)
					continue
				if  int(line.split()[1]) in Remove_bonds:
					continue


				#if int(line.split()[2]) > 45:
				elif int(line.split()[1]) not in Remove_bonds:
				# Any other atoms can just be rewritten as is
					line = line.split()
					new_atoms1 = newAtomNum[oldAtomNum.index(int(line[2]))]
					new_atoms2 = newAtomNum[oldAtomNum.index(int(line[3]))]
					w.write(f'     {curr_bond}      {line[1]}        {new_atoms1}     {new_atoms2}\n')
					curr_bond += 1


		for line in lines[9932:17870]:
			# Remove adsorbate atoms
				print(line.split())
				if len(line.split()) <= 2:
					w.write(line)
					continue
				if  int(line.split()[1]) in Remove_angles:
					continue

				elif int(line.split()[1]) not in Remove_angles:
				# Any other atoms can just be rewritten as is
					line = line.split()
					new_atoms1 = newAtomNum[oldAtomNum.index(int(line[2]))]
					new_atoms2 = newAtomNum[oldAtomNum.index(int(line[3]))]
					new_atoms3 = newAtomNum[oldAtomNum.index(int(line[4]))]

					w.write(f'     {curr_angle}      {line[1]}        {new_atoms1}     {new_atoms2}      {new_atoms3}\n')
					curr_angle += 1

		for line in lines[17870: 25874]:
			# Remove adsorbate atoms
				print(line.split())
				if len(line.split()) <= 2:
					w.write(line)
					continue
				if  int(line.split()[1]) in Remove_dihedrals:
					continue

				elif int(line.split()[1]) not in Remove_dihedrals:
					line = line.split()
					#line1[0] = str(curr_di)
					#line1 = "    ".join(line1)
					#w.write(f'{line1}\n')
					new_atoms1 = newAtomNum[oldAtomNum.index(int(line[2]))]
					new_atoms2 = newAtomNum[oldAtomNum.index(int(line[3]))]
					new_atoms3 = newAtomNum[oldAtomNum.index(int(line[4]))]
					new_atoms4 = newAtomNum[oldAtomNum.index(int(line[5]))]
					w.write(f'     {curr_di}      {line[1]}        {new_atoms1}     {new_atoms2}      {new_atoms3}      {new_atoms4}\n')
					curr_di += 1
		# impropers dont change
		for line in lines[25874:]:
						# Remove adsorbate atoms
			print(line.split())
			if len(line.split()) <= 2:
				w.write(line)
				continue
				#if  int(line.split()[1]) in Remove_dihedrals:
				#	continue

			else:
				line = line.split()
					#line1[0] = str(curr_di)
					#line1 = "    ".join(line1)
					#w.write(f'{line1}\n')
				new_atoms1 = newAtomNum[oldAtomNum.index(int(line[2]))]
				new_atoms2 = newAtomNum[oldAtomNum.index(int(line[3]))]
				new_atoms3 = newAtomNum[oldAtomNum.index(int(line[4]))]
				new_atoms4 = newAtomNum[oldAtomNum.index(int(line[5]))]
				w.write(f'     {line[0]}      {line[1]}        {new_atoms1}     {new_atoms2}      {new_atoms3}      {new_atoms4}\n')
			#w.write(line)

		print(oldAtomNum[453:465])
		print(newAtomNum[453:465:])


		print(oldAtomNum[3600:3616])
		print(newAtomNum[3600:3616])


		print(len(newAtomNum))
		print(len(oldAtomNum))

		new_atoms1index = oldAtomNum.index(int('4348'))
		new_atoms1 = newAtomNum[new_atoms1index]
		print(new_atoms1index)
		print(new_atoms1)