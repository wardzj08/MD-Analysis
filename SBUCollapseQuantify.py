import numpy as np
import pandas as pd

num_atoms = 3648
snp = 50
skips = num_atoms*snp + 9*(snp+1)
input_file = './data/dump.BOYDRELAX.lammpstrj'
atoms = pd.read_csv(input_file, skiprows = skips, nrows = num_atoms, error_bad_lines = False, names = ['Atom','mol', 'type', 'x', 'y', 'z'], delim_whitespace = True)
print(len(atoms))
print(atoms)
       # df = pd.read_csv(input_file, skiprows = skips, nrows = num_atoms, error_bad_lines = False, names = ['Atom','mol', 'type', 'element', 'x', 'y', 'z'], delim_whitespace=True)
mu3O = atoms[atoms['type'] == 4]
print(len(mu3O))
def findMolPairsWithinDistance(ATOM1DF, ATOM2DF, cutoff, oneBondPerAtom1 = False):
        df1_repeated = pd.concat([ATOM1DF[['Atom','x', 'y', 'z']]]*ATOM2DF.shape[0], ignore_index=True)
        df1_repeated.columns = ['Atom1', '1x', '1y', '1z']
        # ATOM2 repeated n times for each atom before the next atom is added to the list
        df2_repeated = ATOM2DF[['Atom', 'x', 'y', 'z']].loc[ATOM2DF.index.repeat(ATOM1DF.shape[0])].reset_index(drop=True)
        df2_repeated.columns = ['Atom2', '2x', '2y', '2z']
        # Store just the coordinates as arrays
        Loc1 = np.array(df1_repeated[['1x', '1y', '1z']])
        Loc2 = np.array(df2_repeated[['2x', '2y', '2z']])

        # Vectorized distance calculation
        # L1,L2 are the 2 arrays of locations for atom1 and atom2 respectivley
        calcDist = lambda L1,L2 : np.sqrt((L1[:, 0] - L2[:, 0]) ** 2 + (L1[:, 1] - L2[:, 1]) ** 2 + (L1[:, 2] - L2[:, 2]) ** 2)
        distances = pd.DataFrame(calcDist(Loc1,Loc2), columns = ['Distance'])

       # Overall data frame of ATOM1-ATOM2
        Overall = pd.concat([df1_repeated, df2_repeated, distances], axis=1).reset_index(drop=True)
        #print('Overall')
        #print(Overall.shape)
        #print(Overall[:5]) # Can check that distances match coordinate values and that atom numbers are coorect by comparing atom numbers in this df to atom numbers in the dump file

        # if oneBondPerAtom1 then finds the minimum distance betweeen atoms for each unique ATOM1 and ignores all other pairings
        # Then checks the remaining pairs (1 for each ATOM1) against the cutoff distance
        # Might be useful for 0 loading step to avoid many acetones being hydrogen bonded to 1 mu3OH (artificially limit number of bonds possible)
        # This would be a pure nearest neighbor implementation
        # Also, with finite loading, setting oneBondPerAtom1 doesn't seem to matter for the result, which makes sense
        if oneBondPerAtom1:
            # This filters out all ATOM2s but the one closest to ATOM1
            atom1NearestNeighbor = Overall.groupby('Atom1')['Distance'].min()
            # Checks against the cutoff value, storing only the pairs that are within that distance
            BondsOccuring = pd.DataFrame({'Distance': atom1NearestNeighbor[atom1NearestNeighbor.values <= cutoff]})
            final = Overall[(Overall['Atom1'].isin(list(BondsOccuring.index))) & (Overall['Distance'].isin(BondsOccuring['Distance']))]

        # Else just check the pair distances against the cutoff distance
        # This has the possibility for more than 1 ATOM2 to be paired with 1 ATOM1 or vice versa
        else:
            #print("final")
            #print(Overall)
            #print(Overall.sort_values('Distance'))
            final = Overall[Overall['Distance'] < cutoff]

        #print('Number of bonded:', final.shape[0])
        return final


#broken_groups =
k = 0
count_broken_groups = 0
for i in range(0,len(mu3O)):
	#print(mu3O.iloc[[i]])
	O = pd.DataFrame(mu3O.iloc[[i]])
	print(O.iloc[0]['Atom'])

	# Check if current atom is already in a "broken group". If it is, it should not be counted again
	if count_broken_groups > 0:
		if O.iloc[0]['Atom'] in list(broken_groups['Atom2']):
			continue
	others = mu3O.drop(mu3O.index[i])
	print(len(others))
	grp = findMolPairsWithinDistance(O, others, 0.5)
	print(grp)
	if len(grp) > 1:
		count_broken_groups += 1
		grp['Group Number'] = k
		if k == 0:
			broken_groups = grp
		else :
			broken_groups = broken_groups.append(grp)
		k += 1

pd.set_option("display.max_columns", None)
print(broken_groups)
print(count_broken_groups)
print(len(broken_groups))

