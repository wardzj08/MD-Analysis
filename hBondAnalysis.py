import pandas as pd
import numpy as np
from itertools import product
pd.set_option('display.max_columns', None)

# variables to be changed based on input
input_file = './data/dump.productionAce1load.lammpstrj'
snapshots_to_read = 1
num_MOF_atoms = 3648 # Atoms in the MOF (UIO66)
MOF_Formula = {'c': 48, 'H': 28, 'O': 32, 'Zr': 6} # Dictionary form of the molecular formula for the MOF
loading = 1 # level of loading  (diffusing molecules per primative cell)
atoms_per_adsorbate = 4 # number of atoms in the diffusing molecule

# Number of atoms = MOFAtoms + AtomsPerAdsorbate * AdsorbateAtomsPerPrimativeCell(loading) * numPrimativeCells(MOFATOMS/ATOMSPERPRIMATIVECELL)
num_atoms = int(num_MOF_atoms + loading * atoms_per_adsorbate * num_MOF_atoms/(sum(MOF_Formula.values())))
print('Atoms in system: ', num_atoms)

MOF_atoms = [444] # MOF molecule number(s) in the LAMMPSTRJ file
mu3_O = 4 # 'Type' of molecule that represents mu3 Oxygens in UiO66
MOF_H = 2 # 'Type' of molecule that represents hydrogens in the MOF (UiO66)
mu3OHBondDist = 1.2 # angstrom distance for the mu3OH bond
hydrogenBondDist = 3 # Acetone-Hydrogen hydrogen bonding distance

# Determines if  atom1s are within cutoff (bonding) distance of atom2s
# for finding mu3 hydrogens, pass dataframe of hydrogens in as atom1 and mu3 oxygens as atom2
# for finding acetone oxygens that are hydrogen bonding to mu3 hydrogens, pass acetone oxygens as atom1 and mu3 hydrogens (the result of the above calc) as atom2
def findMolPairsWithinDistance(ATOM1DF, ATOM2DF, cutoff, oneBondPerAtom1 = False):
        # Create a dataframe of ATOM1-ATOM2 pairings
        # Every ATOM1 paired with every ATOM2
        #pairs = pd.DataFrame(list(product(ATOM1DF['Atom'], ATOM2DF['Atom'])), columns = ['Atom1', 'Atom2'])
     #   print('Number of ATOM1-ATOM2 pairings (should equal number of ATOM1s * number of ATOM2s):', len(pairs))

        # By concatanating the coordinates of ATOM1 k times, a list that matches the number of ATOM1-ATOM2 pairings is
        # created, where k is the number of ATOM2s. This list goes through every ATOM1 atom 1 time before repeating the list again k times

        # In contrast, the ATOM2 coordinate list is repeated n times for each ATOM1 before moving to the next atom
        # where n is the number of ATOM1s in the system

        # By joining these 2 lists with the pairs of atom numbers, a dataframe is formed that has the format:
        # ATOM1 #; ATOM2 #; ATOM1 X coord; ATOM1 Y coord; ATOM1 Z coord; ATOM2 X coord; ATOM2 Y coord; ATOM2 Z coord;
        # This dataframe is also joined by a Distance list, which is the calculated distance between the two atoms atoms in the given row
        # Calculations can be checked easily by printing the Overall dataframe, which has atom#,x,y,z for the two atom types as well as the distance betweeen the two

        # ATOM1 coordinates repeated end on end for the number of ATOM2s there are in the system
        df1_repeated = pd.concat([ATOM1DF[['Atom','x', 'y', 'z']]]*ATOM2DF.shape[0], ignore_index=True)
        df1_repeated.columns = ['Atom1', '1x', '1y', '1z']
        # ATOM2 repeated n times for each atom before the next atom is added to the list.
        df2_repeated = ATOM2DF[['Atom', 'x', 'y', 'z']].loc[ATOM2DF.index.repeat(ATOM1DF.shape[0])].reset_index(drop=True)
        df2_repeated.columns = ['Atom2', '2x', '2y', '2z']
        print(len(df1_repeated), len(df2_repeated))
        # Store just the coordinates as arrays
        Loc1 = np.array(df1_repeated[['1x', '1y', '1z']])
        Loc2 = np.array(df2_repeated[['2x', '2y', '2z']])

        # Vectorized distance calculation
        distances = pd.DataFrame(np.sqrt((Loc1[:, 0] - Loc2[:, 0]) ** 2 + (Loc1[:, 1] - Loc2[:, 1]) ** 2 + (Loc1[:, 2] - Loc2[:, 2]) ** 2), columns = ['Distance'])

       # Overall data frame of ATOM1-ATOM2
        Overall = pd.concat([df1_repeated, df2_repeated, distances], axis=1, sort=False)

        print('Overall')
        print(Overall.shape)
        #print(Overall[:5]) # Can check that distances match coordinate values and that atom numbers are coorect by comparing atom numbers in this df to atom numbers in the dump file

        # if oneBondPerAtom1 then finds the minimum distance betweeen atoms for each ATOM1 and ignores all other pairings
        # Then checks the remaining pairs (1 for each ATOM1) against the cutoff distance
       # print('test')
        #print(Overall[Overall['Distance'] < cutoff])
        if oneBondPerAtom1:
        # This filters out oxygen molecules but the closest one to the hydrogen
            minimumDistPerMover = Overall.groupby('Atom1')['Distance'].min()
            BondsOccuring = pd.DataFrame({'Distance': minimumDistPerMover[minimumDistPerMover.values <= cutoff]})
            final = Overall[(Overall['Atom1'].isin(list(BondsOccuring.index))) & (Overall['Distance'].isin(BondsOccuring['Distance']))]

        # Else just check the pair distances against the cutoff distance
        # This has the possibility for more than 1 ATOM1 to be paired with 1 ATOM2 or vice versa
        else:
            final = Overall[Overall['Distance'] < cutoff]

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
		raise ValueError('The number of H atoms in the MOF read in from the trajectory file do not match the number of H atoms calculated using the number of atoms in the MOF and the molecular formula of the MOF. Check that your trajectory file and the inputed values are correct.')

    # Number of acetone atoms should equal the number of primative cells time the number of adsorbate atoms per cell
	calculatedAceOatoms = num_MOF_atoms/(sum(MOF_Formula.values())) * loading
	if calculatedAceOatoms != readinAceOatoms:
		raise ValueError('The number of Acetone O atoms read in from the trajectory file do not match the number of acetone atoms calculated. Check that your trajectory file and the inputed values are correct.')

# Overall counts for Acetone oxygen - mu3OH hydrogen bonding and total number of Acetone molecules
Ace_mu3HCount = 0
totalAceCount = 0
for snp in range(snapshots_to_read):
        print('Snapshot Number:', snp + 1)
        # Reads in 1 snapshot at a time
        skips = num_atoms*snp + 9*(snp+1)
        df = pd.read_csv(input_file, skiprows = skips, nrows = num_atoms, error_bad_lines = False, names = ['Atom','mol', 'type', 'element', 'x', 'y', 'z'], delim_whitespace=True)

         # Select oxygen molecules in acetone (so all oxygens not in the MOF)
         # As a check the number of atoms in this should equal the number of Acetone molecules inserted
         # So for 2 loading per primative cell and 32 primative cells, there should be 64 Acetone oxygens
        # selects atoms that are Os and not in the MOF
        Ace_Oxygen = df[(~df['mol'].isin(MOF_atoms)) & (df['element'] == 'O')]
        print('Number of Acetone Oxygen Atoms:', len(Ace_Oxygen))

        # Store all MOF hydrogens and mu3 oxygens
        # Each of these have unique type values so they can be selected using only 1 identifier
        MOFHydrogen = df[(df['type'] == MOF_H)]
        MOFOxygen = df[(df['type'] == mu3_O)]
        print('Number of Hydrogens in MOF:', len(MOFHydrogen))
        print('Number of mu3 Oxygens in MOF:', len(MOFOxygen))

        # Only preform the check on the first snapshot/timestep
        # If the first is correct, the remaining should also be
        # Working on making this a bit more general
        if snp == 0:
            num_atoms_check(num_MOF_atoms, loading, len(MOFHydrogen), len(Ace_Oxygen), MOF_Formula)


        # Finds mu3O-H pairs that are within the bonding cutoff distance
        # Atom1 in mu3OH is oxygen, Atom2 is hydrogen; x1,y1,z1 are coords for oxygen; x2,y2,z2 are coords for hydrogen
        mu3OH = findMolPairsWithinDistance(MOFOxygen, MOFHydrogen, mu3OHBondDist, oneBondPerAtom1 = False)
        print('mu3O-H Pairs:')
        print(mu3OH.shape)
        #print(mu3OH.shape)
        # Some dataframe rearanging to format it for passing back into findMolPairsWithinDistance
        # Select only the hydrogens for acetone hydrogen bonding calculations
        mu3H = mu3OH[['Atom2', '2x', '2y', '2z']]
        # Columns need to be 'Atom', 'x', 'y', 'z'
        mu3H.columns = ['Atom', 'x', 'y', 'z']
        print('mu3Hs:')
        print(mu3H[:5])

        # Finds mu3H-AcetoneO pairs that are within the hydrogen bonding cutoff distance
        # Atom1 in is acetone's oxygen, Atom2 is mu3 Hydrogen; x1,y1,z1 are coords for the oxygen; x2,y2,z2 are coords for hydrogen
        aceOmu3H = findMolPairsWithinDistance(Ace_Oxygen, mu3H, hydrogenBondDist, oneBondPerAtom1 = False)
        print('mu3H - Acetone O Hydrogen Bond Pairs:')
        print(aceOmu3H)

        # Add to count of acetone oxygen atoms that are in the hydrogen bonding distance to mu3 hydrogens
        Ace_mu3HCount += len(aceOmu3H)
        # Count the total number of acetone oxygens (should equal loading*primative cells*number of snapshots at the end)
        totalAceCount += len(Ace_Oxygen)


#print('Values across {0} snapshots:'.format(snapshots_to_read))
#print('Number of hydrogen bonding acetone molecules:', Ace_mu3HCount)
#print('Total number of acetone molecules:',totalAceCount)
#print('Fraction of acetone molecules in hydrogen bonds with mu3OHs:', Ace_mu3HCount/totalAceCount)
