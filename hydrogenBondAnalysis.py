import pandas as pd
import numpy as np
from itertools import product

# read in a frame of the trajectory file
input_file = './data/dump.production.lammpstrj'
masses = {'C': 12.01, 'O': 15.999}#, 'H': 1.007}
snapshots_to_read = 501
num_atoms = 3904 # Total atoms in a single snapshot
num_MOF_atoms = 3648
MOF_Formula = {'C': 48, 'H': 28, 'O': 32, 'Zr': 6}
loading = 2
atoms_per_adsorbate = 4 # atoms in the diffusing molecule
MOF_atoms = [444]
mu3_O = 4 # 'Type' of molecule that represents mu3 ohs
MOF_H = 2 # 'Type' of molecule that represents hydrogens in the MOF (UiO66)

def distanceFormula(Loc1, Loc2):
	SUM = (Loc2[0] - Loc1[0])**2 + (Loc2[1] - Loc1[1])**2 + (Loc2[2]-Loc1[2])**2
	return SUM**.5

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
        if snp == 0:
            num_atoms_check(num_MOF_atoms, loading, len(MOFHydrogen), len(Ace_Oxygen), MOF_Formula)


        # Create a dataframe of mu3 O-H pairings
        # Every hydrogen paired with every mu3 oxygen
        # this is used to check the distance from each mu3 oxygen for each hydrogen

        # Pairing of atom numbers for all MOF Hs and mu3 Os
        pairedHmu3O = pd.DataFrame(list(product(MOFHydrogen['Atom'], MOFOxygen['Atom'])), columns = ['HAtom', 'OAtom'])
        print('Number of H-mu3O pairings (should equal number of MOF Hydrogens * number of mu3 Oxygens):', len(pairedHmu3O))

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
        Odf_repeated = pd.concat([MOFOxygen[['x', 'y', 'z']]]*MOFHydrogen.shape[0], ignore_index=True)
        Odf_repeated.columns = ['Ox', 'Oy', 'Oz']
        # Hydrogens repeate n times for each atom before the next atom is added to the list.
        Hdf_repeated = MOFHydrogen[['x', 'y', 'z']].loc[MOFHydrogen.index.repeat(MOFOxygen.shape[0])].reset_index(drop=True)
        Hdf_repeated.columns = ['Hx', 'Hy', 'Hz']

        # Store just the coordinates as arrays
        HLoc = np.array(Hdf_repeated[['Hx', 'Hy', 'Hz']])
        OLoc = np.array(Odf_repeated[['Ox', 'Oy', 'Oz']])

        # Vectorized distance calculation
        distances = pd.DataFrame(np.sqrt((HLoc[:, 0] - OLoc[:, 0]) ** 2 + (HLoc[:, 1] - OLoc[:, 1]) ** 2 + (HLoc[:, 2] - OLoc[:, 2]) ** 2), columns = ['Distance'])

       # Overall data frame of mu3 oxygen - hydrogens
        mu3OH_Overall = pd.concat([pairedHmu3O, Hdf_repeated, Odf_repeated, distances], axis=1, sort=False)

        pd.set_option('display.max_columns', None)
        print(mu3OH_Overall[:5]) # Can check that distances match coordinate values and that atom numbers are coorect by comparing atom numbers in this df to atom numbers in the dump file

        # Find minimum distance for each Hydrogen
        # This filters out oxygen molecules but the closest one to the hydrogen
        mu3OH_minimumPerHydrogen = mu3OH_Overall.groupby('HAtom')['Distance'].min()

        # dataframe of distances of hydrogens in bonding distance to mu3 oxygens
        # Index of this is the atom number for the hydrogen
        # Index #, distance value serve as a key to find the O-H pair in the overall df
        mu3OHBonds = pd.DataFrame({'Distance': mu3OH_minimumPerHydrogen[mu3OH_minimumPerHydrogen < 0.861362]})

        # Store final dataframe of mu3 Hydrogens (in bonding range with mu3 Oxygen)
        # These will be the values that are checked against the acetone oxygen
        # This has the mu3 Oxygen #, mu3 Hydrogen #, coordinates for each and the distance between them
        final_muH = mu3OH_Overall[(mu3OH_Overall['HAtom'].isin(list(mu3OHBonds.index))) & (mu3OH_Overall['Distance'].isin(mu3OHBonds['Distance']))]
        #print(final_muH)
        print('Number of mu3OH bonds:', final_muH.shape[0])

        # Follow the same idea as finding Hs in bonding distance of mu3 Oxygens
        # Will probably write this as a function as its used twice and is fairly long
        paired_Mu3HAceO = pd.DataFrame(list(product(final_muH['HAtom'], Ace_Oxygen['Atom'])), columns = ['mu3HAtom', 'AceOAtom'])

        AceOdf_repeated = pd.concat([Ace_Oxygen[['x', 'y', 'z']]]*final_muH.shape[0], ignore_index=True)
        AceOdf_repeated.columns = ['AceOx', 'AceOy', 'AceOz']
        mu3Hdf_repeated = final_muH[['Hx', 'Hy', 'Hz']].loc[final_muH.index.repeat(Ace_Oxygen.shape[0])].reset_index(drop=True)
        mu3Hdf_repeated.columns = ['3Hx', '3Hy', '3Hz']

        mu3HLoc = np.array(mu3Hdf_repeated[['3Hx', '3Hy', '3Hz']])
        AceOLoc = np.array(AceOdf_repeated[['AceOx', 'AceOy', 'AceOz']])


        Mu3HAceO_distances = pd.DataFrame(np.sqrt((mu3HLoc[:, 0] - AceOLoc[:, 0]) ** 2 + (mu3HLoc[:, 1] - AceOLoc[:, 1]) ** 2 + (mu3HLoc[:, 2] - AceOLoc[:, 2]) ** 2), columns = ['Distance'])

        Ace_Mu3_HbondOverall = pd.concat([paired_Mu3HAceO, mu3Hdf_repeated, AceOdf_repeated, Mu3HAceO_distances], axis=1, sort=False)

        # Find minimum distance for each Hydrogen
        # This filters out oxygen molecules but the closest one to the hydrogen
        mu3HAceO_minimumPerAceOxygen = Ace_Mu3_HbondOverall.groupby('AceOAtom')['Distance'].min()
        print(len(mu3HAceO_minimumPerAceOxygen[mu3HAceO_minimumPerAceOxygen < 3.0]))
        Ace_mu3HCount += len(mu3HAceO_minimumPerAceOxygen[mu3HAceO_minimumPerAceOxygen <3.0])
        totalAceCount += len(Ace_Oxygen)


print(Ace_mu3HCount)
print(totalAceCount)
print(Ace_mu3HCount/totalAceCount)
# Calculate distances between acetone O and mu3 OH hydrogen
# Maybe only take the minimum for each acetone molecule (so the list does not have distances for every acetone to every mu3 OH group)
# Then count number of acetone O atoms within range (ex 3 angstroms)


# do for each frame, calculate overall fraction of acetone molecules engaged in H bonding

# Potential histogram of distances from mu3 group; ie acetone within 3 angstroms 4 angstroms 5 angstroms etc