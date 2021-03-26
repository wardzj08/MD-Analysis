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
input_file = r'D:/LAMMPSTraj/RoggeHbond/dump.NoAceRogge200k.lammpstrj'
input_file = r'D:/LAMMPSTraj/RoggeHbond/dump.BoydNoAceTransition.lammpstrj'
#input_file = "../../../MDResults/dump.EQL1.lammpstrj"
start_snapShot = 0 # Initial snapshot to read
end_snapShot = 20001 # Final snapshot to read
step = 100 # Number of snapshots to skip between read ins

# Number of atoms in the system. Default = 3648 framework atoms + 4 acetone atoms
num_atoms = 3648

# Linkers of interest
# Two linkers set up so atoms pair in the following order - Ca1-C1-C2-Ca2-C3-C4
# When vizualizing UiO-66, the atoms in C_linker1 should be analogous to those in C_linker2 and the two linkers should form 1 transition window
C_linker1 = [2429, 2293, 2288, 2437, 2300, 2281]
C_linker2 = [3313, 3358, 3240, 3302, 3263, 3250]
C_linker3 = [2418, 2318, 2453, 2408, 2309, 2462]

#C_linker1 = [2404, 2449, 2335, 2398, 2443, 2325]
#C_linker2 = [2434, 2298, 2471, 2427, 2307, 2286]

# Rogge
#C_linker1 = [3549, 3215, 3419, 3551, 3224, 3404]
#C_linker2 = [827, 611, 515, 844, 641, 545]
#C_linker3 = [1765, 1499, 1631, 1743, 1478, 1634]


# Distance function. Though this is a periodic system, no need to apply min. image conventions as this is only looking at atoms that are in the center of the sim. box
calcDist = lambda L1,L2 : np.sqrt((L1[:, 0] - L2[:, 0]) ** 2 + (L1[:, 1] - L2[:, 1]) ** 2 + (L1[:, 2] - L2[:, 2]) ** 2)

# Create a list of the snapshots to use
snapsToRead = range(start_snapShot, end_snapShot, step)
totalSnaps = len(snapsToRead)
print(totalSnaps)
# final distances array will have totalSnaps rows (i.e. number of snapshots used) and len(C_linker1) columns (number of atom pairs to calculate distances from)
# Each row of the array will be filled as snapshots are looped over with the distances for each of the atom pairs
distances12 = np.zeros(shape = (totalSnaps, 2))
distances13 = np.zeros(shape = (totalSnaps, 2))
distances23 = np.zeros(shape = (totalSnaps, 2))

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
        df1 = df[df['Atom'].isin(C_linker1)].set_index('Atom').reindex(C_linker1)
        df2 = df[df['Atom'].isin(C_linker2)].set_index('Atom').reindex(C_linker2)
        df3 = df[df['Atom'].isin(C_linker3)].set_index('Atom').reindex(C_linker3)
        Loc1 = np.array(df1[['x', 'y', 'z']])
        Loc2 = np.array(df2[['x', 'y', 'z']])
        Loc3 = np.array(df3[['x', 'y', 'z']])

        # Averaging inside and outside carbons into 1 location, omit for atom by atom comparison
        L1Avgs = np.zeros(shape = (2, 3))
        L2Avgs = np.zeros(shape = (2, 3))
        L3Avgs = np.zeros(shape = (2, 3))
        for link, out in zip([Loc1, Loc2, Loc3], [L1Avgs, L2Avgs, L3Avgs]):
        	xin = (link[1, 0] + link[2, 0]) * 0.5
        	yin = (link[1, 1] + link[2, 1]) * 0.5
        	zin = (link[1, 2] + link[2, 2]) * 0.5
        	out[0, :] = [xin, yin, zin]

        	xout = (link[4, 0] + link[5, 0]) * 0.5
        	yout = (link[4, 1] + link[5, 1]) * 0.5
        	zout = (link[4, 2] + link[5, 2]) * 0.5
        	out[1, :] = [xout, yout, zout]

        print(Loc1)
        print(L1Avgs)
        print(L2Avgs)

        # Make the distance calculations. Add to the overall array
        #distance = np.array(calcDist(Loc1,Loc2))
        distance12 = np.array(calcDist(L1Avgs,L2Avgs)) # Distances between the average of the two sides
        distance13 = np.array(calcDist(L1Avgs,L3Avgs)) # Distances between the average of the two sides
        distance23 = np.array(calcDist(L2Avgs,L3Avgs)) # Distances between the average of the two sides
        distances12[curr, :] = distance12
        distances13[curr, :] = distance13
        distances23[curr, :] = distance23
        curr = curr + 1
        #print(distances)


# Initial plotting
import matplotlib.pyplot as plt
# Save the array to a csv file. Easier plotting without redoing calculations in file: plotting_WindowSize.py
np.savetxt("BoydNoAceDistancesAveragedInOut12.csv", distances12, delimiter=",")
np.savetxt("BoydNoAceStartDistancesAveragedInOut13.csv", distances13, delimiter=",")
np.savetxt("BoydNoAceStartDistancesAveragedInOut23.csv", distances23, delimiter=",")
plt.plot(range(totalSnaps), distances12[:, 0], label = '1-2Interior')
plt.plot(range(totalSnaps), distances12[:, 1], label = '1-2Exterior')
plt.plot(range(totalSnaps), distances13[:, 0], label = '1-3Interior')
plt.plot(range(totalSnaps), distances13[:, 1], label = '1-3Exterior')
plt.plot(range(totalSnaps), distances23[:, 0], label = '2-3Interior')
plt.plot(range(totalSnaps), distances23[:, 1], label = '2-3Exterior')
plt.legend()

plt.savefig('RoggeNoAceAveragedLinkers.PNG')
