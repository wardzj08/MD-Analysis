# Given a LAMMPS trajectory file of an adsorbate in UiO-66 system, calculates the distances between analagous atoms on linkers forming 1 transtion window. Goal is to
# understand how the transition state changes as an adsorbate jumps between cages. This will explain why rigid models for UiO-66 do not work and additionally why a rigid
# framework could not be designed to properly model this system.

import numpy as np
import pandas as pd


input_file = "../../../MDResults/dump.FULLRC3.lammpstrj"
#input_file = r'D:/LAMMPSTraj/RoggeHbond/dump.RoggeHbondTetraOctaTransition.lammpstrj'
#input_file = r'D:/LAMMPSTraj/RoggeHbond/dump.RoggeHbond25K.lammpstrj'
#input_file = r'D:/LAMMPSTraj/dump.RoggeFullRC3.lammpstrj'
#input_file = "../../../MDResults/dump.FULLRC.lammpstrj" # trajectory file to read from
# No ace in window
#input_file = r'D:/LAMMPSTraj/RoggeHbond/dump.NoAceRogge200k.lammpstrj'
#input_file = r'D:/LAMMPSTraj/RoggeHbond/dump.BoydNoAceTransition.lammpstrj'
#input_file = "../../../MDResults/dump.EQL1.lammpstrj"

start_snapShot = 0 # Initial snapshot to read
end_snapShot = 20001 # Final snapshot to read
step = 100TETRACOMDEN = np.sum(dfTetra['mass']) # Number of snapshots to skip between read ins

# Number of atoms in the system. Default = 3648 framework atoms + 4 acetone atoms
num_atoms = 3652

# Linkers of interest
# Two linkers set up so atoms pair in the following order - Ca1-C1-C2-Ca2-C3-C4
# When vizualizing UiO-66, the atoms in C_linker1 should be analogous to those in C_linker2 and the two linkers should form 1 transition window

FindCOM = [2714, 3633, 3646] # zrs
AceC = 3650


# Distance function. Though this is a periodic system, no need to apply min. image conventions as this is only looking at atoms that are in the center of the sim. box
calcDist = lambda L1,L2 : np.sqrt((L1[:, 0] - L2[:, 0]) ** 2 + (L1[:, 1] - L2[:, 1]) ** 2 + (L1[:, 2] - L2[:, 2]) ** 2)

# Create a list of the snapshots to use
snapsToRead = range(start_snapShot, end_snapShot, step)
totalSnaps = len(snapsToRead)
print(totalSnaps)
# final distances array will have totalSnaps rows (i.e. number of snapshots used) and len(C_linker1) columns (number of atom pairs to calculate distances from)
# Each row of the array will be filled as snapshots are looped over with the distances for each of the atom pairs
distancesAceCen = np.zeros(shape = (totalSnaps, 1))

curr = 0 # tracker for current row
for snp in snapsToRead:
        print('Snapshot Number:', snp + 1)
        # Reads in 1 snapshot at a time
        skips = num_atoms*snp + 9*(snp+1)
        df = pd.read_csv(input_file, skiprows = skips, nrows = num_atoms, error_bad_lines = False, names = ['Atom','mol', 'type', 'x', 'y', 'z'], delim_whitespace=True)

        # Create subsets
        # 1 is the carbon atoms on linker 1
        # 2 is the carbon atoms on linker 2
        # .set_index('Atom').reindex(C_linkerx) is important to prevent reindexing of atoms based on numerical order
        # If this is not used, the two dfs will not be in the needed order to compare atoms
        df1 = df[df['Atom'].isin(FindCOM)]
       	AceAtom = df[df['Atom'] == 3650]
        CALCCOMDIR = lambda pos, mass, sumMass : np.sum(pos * mass) / sumMass

        posCOMWin = []
        DEN = 91.222 * 3
        # Calculate the COM in each dimension for each cage
        for pos in ['x', 'y', 'z']:
            posCOMWin = posCOMWin + [CALCCOMDIR(df1[pos], [91.222, 91.222, 91.222], DEN)]

        #df2 = df[df['Atom'].isin(C_linker2)].set_index('Atom').reindex(C_linker2)
        #df3 = df[df['Atom'].isin(C_linker3)].set_index('Atom').reindex(C_linker3)
        Loc1 = np.array([posCOMWin])
        print(Loc1)
        Loc2 = np.array(AceAtom[['x', 'y', 'z']])
        #Loc3 = np.array(df3[['x', 'y', 'z']])

        print(Loc2)
        # Make the distance calculations. Add to the overall array
        #distance = np.array(calcDist(Loc1,Loc2))
        distanceAC = np.array(calcDist(Loc1,Loc2)) # Distances between the average of the two sides

        distancesAceCen[curr, :] = distanceAC

        curr = curr + 1
        #print(distances)

print(distancesAceCen)
# Initial plotting
import matplotlib.pyplot as plt
# Save the array to a csv file. Easier plotting without redoing calculations in file: plotting_WindowSize.py
np.savetxt("AceCenDist.csv", distancesAceCen, delimiter=",")

plt.plot(range(totalSnaps), distancesAceCen)
plt.legend()

plt.savefig('AceToCenterDist.PNG')
