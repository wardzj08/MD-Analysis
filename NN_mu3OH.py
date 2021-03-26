import numpy as np
import pandas as pd

# Change num_atoms based on acetone in system
num_atoms = 3648
snp = 0
skips = num_atoms*snp + 9*(snp+1)
#

input_file = './data/dump.LARGEAngleTrappEH.lammpstrj'
#input_file = './data/dump.EvenLargerAngle14K.lammpstrj'

# Remove element depending on input type style
atoms = pd.read_csv(input_file, skiprows = skips, nrows = num_atoms, error_bad_lines = False, names = ['Atom','mol', 'type', 'x', 'y', 'z'], delim_whitespace = True)
print(len(atoms))
#print(atoms)
       # df = pd.read_csv(input_file, skiprows = skips, nrows = num_atoms, error_bad_lines = False, names = ['Atom','mol', 'type', 'element', 'x', 'y', 'z'], delim_whitespace=True)
mu3H = atoms[atoms['type'] == 6]
print(mu3H.head())
print(len(mu3H))
Nonmu3OAtoms = atoms
print(Nonmu3OAtoms.head())
print(len(Nonmu3OAtoms))
#Zr = atoms[atoms['type'] == 5]

def findMolPairsWithinDistance(ATOM1DF, ATOM2DF, cutoff, oneBondPerAtom1 = False):
        df1_repeated = pd.concat([ATOM1DF[['Atom', 'type', 'x', 'y', 'z']]]*ATOM2DF.shape[0], ignore_index=True)
        df1_repeated.columns = ['Atom1', 'type1', '1x', '1y', '1z']
        # ATOM2 repeated n times for each atom before the next atom is added to the list
        df2_repeated = ATOM2DF[['Atom', 'type', 'x', 'y', 'z']].loc[ATOM2DF.index.repeat(ATOM1DF.shape[0])].reset_index(drop=True)
        df2_repeated.columns = ['Atom2', 'type2', '2x', '2y', '2z']
        # Store just the coordinates as arrays
        Loc1 = np.array(df1_repeated[['1x', '1y', '1z']])
        Loc2 = np.array(df2_repeated[['2x', '2y', '2z']])
        print(len(Loc1), len(Loc2))
        print(Loc2[:5])
        # Vectorized distance calculation
        # L1,L2 are the 2 arrays of locations for atom1 and atom2 respectivley
        calcDist = lambda L1,L2 : np.sqrt((L1[:, 0] - L2[:, 0]) ** 2 + (L1[:, 1] - L2[:, 1]) ** 2 + (L1[:, 2] - L2[:, 2]) ** 2)
        distances = pd.DataFrame(calcDist(Loc1,Loc2), columns = ['Distance'])
        print(len(distances))
        #bounds = np.array([41.401, 41.401, 41.401])
        #min_dists = np.min(np.dstack(((Loc1 - Loc2) % bounds, (Loc2 - Loc1) % bounds)), axis = 2)
        #distances = np.sqrt(np.sum(min_dists ** 2, axis = 1))
        #print(distances)

# min image convention distance funtion
        def distance(x0, x1, dimensions):
            delta = np.abs(x0 - x1)
            delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
            return np.sqrt((delta ** 2).sum(axis=-1))

        distances = pd.DataFrame(distance(Loc1, Loc2, np.array([41.4008, 41.4008, 41.4008])), columns = ['Distance'])
       # print(distances1)
       # print(distances)

               # Overall data frame of ATOM1-ATOM2
        Overall = pd.concat([df1_repeated, df2_repeated, distances], axis=1).reset_index(drop=True)
        #print('Overall')
        print(Overall.head())
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
            Overall = Overall.sort_values('Distance')
            final = Overall[Overall['Distance'] < cutoff]

        print('Number of bonded:', final.shape[0])
        return final





grp = findMolPairsWithinDistance(mu3H, Nonmu3OAtoms, 2.5)

		#grp['SBUNUM'] = k
		#O['SBUNUM'] = k
		#brk1 = brk1.append({'Atom': O.iloc[0]['Atom'], 'SBUNUM': k}, ignore_index = True)
		#brk1 = brk1.append({'Atom': grp['Atom2'], 'SBUNUM': [k]*len(grp)}, ignore_index = True)
		#if k == 0:
		#	brk1 = pd.DataFrame(O[['Atom', 'SBUNUM']])
		#	brk1 = brk1.append(grp[['Atom2', 'SBUNUM']], ignore_index = True)
		#	broken_groups = grp
		#else :
		#	broken_groups = broken_groups.append(grp)
		#	brk1 = brk1.append(O[['Atom', 'SBUNUM']], ignore_index = True)
		#	brk1 = brk1.append(grp[['Atom2', 'SBUNUM']], ignore_index = True)
		#k += 1

pd.set_option("display.max_columns", None)
print(grp.tail(50))
print(len(grp))
print(grp['type2'].value_counts())
print(grp.groupby('type2', as_index=False)['Distance'].mean())
#print(grp[grp['type2'] == 1])
#print(broken_groups)
#print(count_broken_groups)
#print(len(broken_groups))
#print(brk1)
