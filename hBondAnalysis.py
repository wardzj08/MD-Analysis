import pandas as pd
import numpy as np
# Shows all columns of dataframe instead of the front and end only
pd.set_option('display.max_columns', None)

# variables to be changed based on input
input_file = './data/dump.productionAce2load.lammpstrj' # file to read trajectory from
snapshots_to_read = 501 # Number of timesteps to read from traj file
num_MOF_atoms = 3648 # Atoms in the MOF (UIO66)
MOF_Formula = {'c': 48, 'H': 28, 'O': 32, 'Zr': 6} # Dictionary form of the molecular formula for the MOF
loading = 2 # level of loading  (diffusing molecules per primative cell)
atoms_per_adsorbate = 4 # number of atoms in the diffusing molecule

# total atoms in system = MOFAtoms + AtomsPerAdsorbate * AdsorbateAtomsPerPrimativeCell(loading) * numPrimativeCells(MOFATOMS/ATOMSPERPRIMATIVECELL)
num_atoms = int(num_MOF_atoms + loading * atoms_per_adsorbate * num_MOF_atoms / (sum(MOF_Formula.values())))
print('Atoms in system: ', num_atoms)

MOF_atoms = [444] # MOF molecule number(s) in the LAMMPSTRJ file
mu3_O = 4 # 'Type' of atom that represents mu3 Oxygens in UiO66
MOF_H = 2 # 'Type' of atom that represents hydrogens in the MOF (UiO66)
Ace_O = 6 # 'Type' of acetone oxygen atom
mu3OHBondDist = 1.2 # angstrom distance for the mu3OH bond
hydrogenBondDist = 3.3 # Acetone-Hydrogen hydrogen bonding distance

# Determines if  atom1s are within cutoff (bonding) distance of atom2s
# for finding mu3 hydrogens, pass dataframe of hydrogens in as atom1 and mu3 oxygens as atom2
# for finding acetone oxygens that are hydrogen bonding to mu3 hydrogens, pass acetone oxygens as atom1 and mu3 hydrogens (the result of the above calc) as atom2
# Don't think the order passed in should have an affect
def findMolPairsWithinDistance(ATOM1DF, ATOM2DF, cutoff, oneBondPerAtom1 = False):

        # Create a dataframe of ATOM1-ATOM2 pairings
        # Every ATOM1 paired with every ATOM2

        # By concatanating the coordinates of ATOM1 k times, a list that matches the number of ATOM1-ATOM2 pairings is
        # created, where k is the number of ATOM2s. This list goes through every ATOM1 atom 1 time before repeating the list again k times

        # In contrast, each individual ATOM2 is repeated n times before moving to the next atom
        # where n is the number of ATOM1s in the system

        # By joining these 2 dataframes a dataframe is formed that has the format:
        # ATOM1#; ATOM1 X coord; ATOM1 Y coord; ATOM1 Z coord; ATOM2#; ATOM2 X coord; ATOM2 Y coord; ATOM2 Z coord;
        # This dataframe is later joined by a cooresponding distance vector
        # Calculations can be checked easily by printing the Overall dataframe
        # The atom numbers match those in the trajectory file, so can be cross-refrenced there

        # ATOM1 coordinates repeated end on end for the number of ATOM2s there are in the system
        df1_repeated = pd.concat([ATOM1DF[['Atom','x', 'y', 'z']]]*ATOM2DF.shape[0], ignore_index=True)
        df1_repeated.columns = ['Atom1', '1x', '1y', '1z']
        # ATOM2 repeated n times for each atom before the next atom is added to the list
        df2_repeated = ATOM2DF[['Atom', 'x', 'y', 'z']].loc[ATOM2DF.index.repeat(ATOM1DF.shape[0])].reset_index(drop=True)
        df2_repeated.columns = ['Atom2', '2x', '2y', '2z']
 #       print(len(df1_repeated), len(df2_repeated))
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

# Function that checks that the number of certain atoms match expectations
# Checks the number of hydrogen atoms in the MOF against the expected value
# Checks number of acetone molecules against the expected value at the given loading
def num_atoms_check(num_MOF_atoms, loading, readinHatoms, readinAceOatoms, MOF_Formula):

	calculatedHatoms = MOF_Formula['H'] * num_MOF_atoms/(sum(MOF_Formula.values()))
	# Check if the number of Hydrogen atoms in the MOF read in mataches the calculated value
	if calculatedHatoms != readinHatoms:
		raise ValueError('The number of H atoms in the MOF read in from the trajectory file do not match the number of H atoms calculated using the number of atoms in the MOF and the molecular formula of the MOF. Check input.')

    # Number of acetone atoms should equal the number of primative cells time the number of adsorbate atoms per cell
	calculatedAceOatoms = num_MOF_atoms/(sum(MOF_Formula.values())) * loading
	if calculatedAceOatoms != readinAceOatoms:
		raise ValueError('The number of Acetone oxygen atoms read in from the trajectory file do not match the number of acetone atoms calculated. Check that your trajectory file and the inputed values are correct.')

Ace_mu3HCount = 0 # Count of hydrogen bonds occuring between mu3H and Acetone
totalAceCount = 0 # Count of total acetone molecules (also equals loading*primativeCells*numSnapshots)
totalmu3HCount = 0 # Count of total mu3H molecules
run_fractionsAce = [] # list of individual snapshots hydrogen bonding fractions in relation to acetone molecules
run_fractionsmu3H = [] # list of individual snapshots hydrogen bonding fractions in relation to mu3 hydrogens
for snp in range(snapshots_to_read):
        print('\nSnapshot Number:', snp + 1)
        # Reads in 1 snapshot at a time
        skips = num_atoms*snp + 9*(snp+1)
        df = pd.read_csv(input_file, skiprows = skips, nrows = num_atoms, error_bad_lines = False, names = ['Atom','mol', 'type', 'element', 'x', 'y', 'z'], delim_whitespace=True)

         # Select oxygen molecules in acetone (so all oxygens not in the MOF)
         # As a check the number of atoms in this should equal the number of Acetone molecules inserted
         # So for 2 loading per primative cell and 32 primative cells, there should be 64 Acetone oxygens
        # selects atoms that are Os and not in the MOF
        Ace_Oxygen = df[df['type'] == Ace_O]
        print('Number of Acetone Oxygen Atoms:', len(Ace_Oxygen))

        # Store all MOF hydrogens and mu3 oxygens
        # Each of these have unique type values so they can be selected using only 1 identifier
        MOFHydrogen = df[(df['type'] == MOF_H)]
        MOFOxygen = df[(df['type'] == mu3_O)]
        print('Number of Hydrogen in MOF:', len(MOFHydrogen))
        print('Number of mu3-Oxygen in MOF:', len(MOFOxygen))

        # Only preform the check on the first snapshot/timestep
        # If the first is correct, the remaining should also be
        # Working on making this a bit more general
        if snp == 0:
            num_atoms_check(num_MOF_atoms, loading, len(MOFHydrogen), len(Ace_Oxygen), MOF_Formula)


        # Finds mu3O-H pairs that are within the bonding cutoff distance
        # Atom1 in mu3OH is oxygen, Atom2 is hydrogen; x1,y1,z1 are coords for oxygen; x2,y2,z2 are coords for hydrogen
        mu3OH = findMolPairsWithinDistance(MOFOxygen, MOFHydrogen, mu3OHBondDist, oneBondPerAtom1 = False)
        #print('mu3O-H Pairs:')
#        print(mu3OH.shape)
        #print(mu3OH.shape)

        # Some dataframe rearanging to format it for passing back into findMolPairsWithinDistance
        # Select only the hydrogens for acetone hydrogen bonding calculations
        mu3H = mu3OH[['Atom2', '2x', '2y', '2z']]
        # Columns need to be 'Atom', 'x', 'y', 'z'
        mu3H.columns = ['Atom', 'x', 'y', 'z']

        # show dataframe of just mu3 hydrogen
        print('Number of mu3-Hs: ', len(mu3H))
       # print('mu3H df: ')
       #print(mu3H[:5])

        # Finds mu3H-AcetoneO pairs that are within the hydrogen bonding cutoff distance
        # Atom1 in is acetone's oxygen, Atom2 is mu3 Hydrogen; x1,y1,z1 are coords for the oxygen; x2,y2,z2 are coords for hydrogen
        aceOmu3H = findMolPairsWithinDistance(Ace_Oxygen, mu3H, hydrogenBondDist, oneBondPerAtom1 = False)

        # show the number of hydrogen bonding pairs
        print('Number of mu3OH - Acetone Hydrogen Bond Pairs: ', len(aceOmu3H))
        # Show pairs of mu3 hydrogen and acetone oxygen that are in hydrogen bonding distance
        #print('Acetone-mu3OH Bond df: ')
        #print(aceOmu3H)

        # Add to count of acetone oxygen atoms that are in the hydrogen bonding distance to mu3 hydrogens
        run_fractionsAce = run_fractionsAce + [len(aceOmu3H)/len(Ace_Oxygen)]
        Ace_mu3HCount += len(aceOmu3H)
        # Count the total number of acetone oxygens (should equal loading*primative cells*number of snapshots at the end)
        totalAceCount += len(Ace_Oxygen)

        # add to fraction of mu3H hydrogen bonding list
        run_fractionsmu3H = run_fractionsmu3H + [len(aceOmu3H)/len(mu3H)]
        # Count the total number of mu3Hs in the system (128 each snapshot)
        totalmu3HCount += len(mu3H)


print('\nValues across {0} snapshots:'.format(snapshots_to_read))
print('Number of hydrogen bonds molecules:', Ace_mu3HCount)
print('Total number of acetone molecules:', totalAceCount)
print('Fraction of acetone molecules in hydrogen bonds with mu3OHs:', round(float(Ace_mu3HCount)/float(totalAceCount), 6))
print(f'Acetone average {round(np.mean(run_fractionsAce), 6)} + std deviation: {round(np.std(run_fractionsAce), 6)}')

print('\nTotal number of mu3H molecules:', totalmu3HCount)
print('Fraction of mu3H molecules in hydrogen bonds with acetone:', round(float(Ace_mu3HCount)/float(totalmu3HCount), 6))
print(f'mu3H average {round(np.mean(run_fractionsmu3H), 6)} + std deviation: {round(np.std(run_fractionsmu3H), 6)}')