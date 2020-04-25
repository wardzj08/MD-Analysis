import pandas as pd
import numpy as np
from itertools import product

# read in a frame of the trajectory file
input_file = 'dump.production.lammpstrj'
masses = {'C': 12.01, 'O': 15.999}#, 'H': 1.007}
snapshots_to_read = 501
num_atoms = 3904 # Total atoms in a single snapshot
num_MOF_atoms = 3648
MOF_Formula = {'c': 48, 'H': 28, 'O': 32, 'Zr': 6}
loading = 2
atoms_per_adsorbate = 4 # atoms in the diffusing molecule
MOF_atoms = [444]
mu3_O = 4 # 'Type' of molecule that represents mu3 ohs
MOF_H = 2 # 'Type' of molecule that represents hydrogens in the MOF (UiO66)

# Determines if molecule 1 atoms are within cutoff (bonding) distance of atom2
# for finding mu3 hydrogens, pass dataframe of hydrogens in as atom1 and mu3 oxygens as atom2
# for finding acetone oxygens that are hydrogen bonding to mu3 hydrogens, pass acetone oxygens as atom1 and mu3 hydrogens (the result of the above calc) as atom2
# I think passing the atom with a lower total count as atom1 and the one with more atoms as atom2
# Basically, the one with less is the limiting atom in the distance calculation. We only want 1 atom1-atom2 pairing in the end, so passing the lesser count as 1
# results in sorting by atom1s #, so there are no issues with having the atom appear as the closest to more than 1 atom2s
def findMolPairsWithinDistance(ATOM1DF, ATOM2DF, cutoff):
        # Create a dataframe of mu3 O-H pairings
        # Every hydrogen paired with every mu3 oxygen
        # this is used to check the distance from each mu3 oxygen for each hydrogen

        # Pairing of atom numbers for all MOF Hs and mu3 Os
        pairs = pd.DataFrame(list(product(ATOM1DF['Atom'], ATOM2DF['Atom'])), columns = ['Atom1', 'Atom2'])
     #   print('Number of H-mu3O pairings (should equal number of MOF Hydrogens * number of mu3 Oxygens):', len(pairedHmu3O))

        # By concatanating the coordinates of the Mu3 Oxygen k times, a list that matches the number of mu3O-H pairings is
        # created, where k is the number of hydrogen atoms in the MOF. This list goes through every mu3 Oxygen atom 1 time
        # before repeating the list again k times

        # In contrast, the MOF hydrogen coordinate list is repeated n times for each hydrogen atom before moving to the next atom
        # where n is the number of mu3 Oxygens in the system

        # By joining these 2 lists with the pairs of atom numbers, a dataframe is formed that has the format:
        # Hydrogen Atom #; mu3 Oxygen Atom #; Hydrogen X coord; Hydrogen Y coord; Hydrogen Z coord; Oxygen X coord; Oxygen Y coord; Oxygen Z coord;
        # This dataframe is also joined by a Distance series, which is the calculated distance between the oxygen and hydrogen atoms in the given row
        # Calculations can be checked easily by printing the mu3OH_Overall dataframe, which has atom#,x,y,z for the hydrogen and oxygen atom as well as the distance betweeen the two

        # mu3 oxygen coordinates repeated end on end for the number of hydrogens there are (should be 896)
        df1_repeated = pd.concat([ATOM1DF[['x', 'y', 'z']]]*ATOM2DF.shape[0], ignore_index=True)
        df1_repeated.columns = ['1x', '1y', '1z']
        # Hydrogens repeate n times for each atom before the next atom is added to the list.
        df2_repeated = ATOM2DF[['x', 'y', 'z']].loc[ATOM2DF.index.repeat(ATOM1DF.shape[0])].reset_index(drop=True)
        df2_repeated.columns = ['2x', '2y', '2z']
        print(len(df1_repeated), len(df2_repeated))
        # Store just the coordinates as arrays
        Loc1 = np.array(df1_repeated[['1x', '1y', '1z']])
        Loc2 = np.array(df2_repeated[['2x', '2y', '2z']])

        # Vectorized distance calculation
        distances = pd.DataFrame(np.sqrt((Loc1[:, 0] - Loc2[:, 0]) ** 2 + (Loc1[:, 1] - Loc2[:, 1]) ** 2 + (Loc1[:, 2] - Loc2[:, 2]) ** 2), columns = ['Distance'])

       # Overall data frame of mu3 oxygen - hydrogens
        Overall = pd.concat([pairs, df1_repeated, df2_repeated, distances], axis=1, sort=False)

        pd.set_option('display.max_columns', None)
        #print(Overall[:5]) # Can check that distances match coordinate values and that atom numbers are coorect by comparing atom numbers in this df to atom numbers in the dump file

        # Find minimum distance for each Hydrogen
        # This filters out oxygen molecules but the closest one to the hydrogen
        minimumDistPerMover = Overall.groupby('Atom1')['Distance'].min()

        # dataframe of distances of hydrogens in bonding distance to mu3 oxygens
        # Index of this is the atom number for the hydrogen
        # Index #, distance value serve as a key to find the O-H pair in the overall df
        BondsOccuring = pd.DataFrame({'Distance': minimumDistPerMover[minimumDistPerMover < cutoff]})

        # Store final dataframe of mu3 Hydrogens (in bonding range with mu3 Oxygen)
        # These will be the values that are checked against the acetone oxygen
        # This has the mu3 Oxygen #, mu3 Hydrogen #, coordinates for each and the distance between them
        final = Overall[(Overall['Atom1'].isin(list(BondsOccuring.index))) & (Overall['Distance'].isin(BondsOccuring['Distance']))]
        #print(final_muH)
        #print('Number of bonded:', final.shape[0])
        return final
# Function that checks that the number of certain atoms match expectations
# Checks the number of hydrogen atoms in the MOF against the expected value
# Checks number of acetone molecules against the expected value
def num_atoms_check(num_MOF_atoms, loading, readinHatoms, readinAceOatoms, MOF_Formula):

	calculatedHatoms = MOF_Formula['H'] * num_MOF_atoms/(sum(MOF_Formula.values()))
	# Check if the number of Hydrogen atoms in the MOF read in mataches the calculated value
	if calculatedHatoms != readinHatoms:
		raise ValueError('The number of hydrogen atoms in the MOF read in from the trajectory file do not match the number of hydrogen atoms calculated using the number of atoms in the MOF and the molecular formula of the MOF. Check that your trajectory file and the inputed values are correct.')

    # Number of acetone atoms should equal the number of primative cells time the number of adsorbate atoms per cell
	calculatedAceOatoms = num_MOF_atoms/(sum(MOF_Formula.values())) * loading
	if calculatedAceOatoms != readinAceOatoms:
		raise ValueError('The number of oxygen atoms in acetone read in from the trajectory file do not match the number of acetone atoms calculated using the . Check that your trajectory file and the inputed values are correct.')

# Overall counts for Acetone oxygen - mu3OH hydrogen bonding and total number of Acetone molecules
Ace_mu3HCount = 0
totalAceCount = 0
for snp in range(snapshots_to_read):
        print('Snapshot Number:', snp + 1)
        # Reads in 1 snapshot at a time
        if snp == 0:
           df = pd.read_csv(input_file, skiprows = 9, nrows = num_atoms, names = ['Atom', 'mol', 'type', 'element', 'x', 'y', 'z'], usecols=[0,1,2,3,4,5,6], delim_whitespace=True)
        else:
            skips = num_atoms*snp + 9*(snp+1)
            df = pd.read_csv(input_file, skiprows = skips, nrows = num_atoms, error_bad_lines = False, names = ['Atom','mol', 'type', 'element', 'x', 'y', 'z'], usecols=[0,1,2,3,4,5,6], delim_whitespace=True)

         # Select oxygen molecules in acetone (so all oxygens not in the MOF)
         # something like if O and if mol # != mol number MOF
         # As a check the number of atoms in this should equal the number of Acetone molecules inserted
         # So for 2 loading per primative cell and 32 primative cells, there should be 64 Acetone oxygens
        Ace_Oxygen = df[(~df['mol'].isin(MOF_atoms)) & (df['element'] == 'O')]
        print('Number of Acetone Oxygen Atoms:', len(Ace_Oxygen))

        # Store all mof hydrogens and mu oxygens
        MOFHydrogen = df[df['type'] == MOF_H]
        MOFOxygen = df[(df['type'] == mu3_O)]
        print('Number of Hydrogens in MOF:', len(MOFHydrogen))
        print('Number of mu3 Oxygens in MOF:', len(MOFOxygen))

        # Only preform the check on the first snapshot/timestep
        # If the first is correct, the remaining should also be
        # Working on making this a bit more general
        if snp == 0:
            num_atoms_check(num_MOF_atoms, loading, len(MOFHydrogen), len(Ace_Oxygen), MOF_Formula)


        # Finds mu3O-H pairs that are within the bonding cutoff distance
        # Atom1 in mu3OH is hydrogen, Atom2 is Oxygen; x1,y1,z1 are coords for hydrogen; x2,y2,z2 are coords for oxygen
        mu3OH = findMolPairsWithinDistance(MOFHydrogen, MOFOxygen, 1.5)
        print(mu3OH.shape)
        # Some dataframe rearanging to format it for passing back into findMolPairsWithinDistance
        mu3H = mu3OH[['Atom1', '1x', '1y', '1z']]
        #print(mu3H.shape)
        # Only need the hydrogen atoms for further calculations. Columns need to be 'Atom', 'x', 'y', 'z'
        mu3H.columns = ['Atom', 'x', 'y', 'z']

        # Finds mu3H-AcetoneO pairs that are within the hydrogen bonding cutoff distance
        # Atom1 in aceOmu3H is acetone's oxygen, Atom2 is mu3Hydrogens; x1,y1,z1 are coords for the oxygen; x2,y2,z2 are coords for hydrogen
        aceOmu3H = findMolPairsWithinDistance(Ace_Oxygen, mu3H, 3.5)

        # Add to count of acetone oxygen atoms that are in the hydrogen bonding distance to mu3 hydrogens
        Ace_mu3HCount += len(aceOmu3H)
        # Count the total number of acetone oxygens (should equal loading*primative cells*number of snapshots at the end)
        totalAceCount += len(Ace_Oxygen)


print('Values across {0} snapshots:', {snapshots_to_read})
print('Number of hydrogen bonding acetone molecules:', Ace_mu3HCount)
print('Total number of acetone molecules:',totalAceCount)
print('Fraction of acetone molecules in hydrogen bonds with mu3OHs:', Ace_mu3HCount/totalAceCount)
# Calculate distances between acetone O and mu3 OH hydrogen
# Maybe only take the minimum for each acetone molecule (so the list does not have distances for every acetone to every mu3 OH group)
# Then count number of acetone O atoms within range (ex 3 angstroms)


# do for each frame, calculate overall fraction of acetone molecules engaged in H bonding

# Potential histogram of distances from mu3 group; ie acetone within 3 angstroms 4 angstroms 5 angstroms etc