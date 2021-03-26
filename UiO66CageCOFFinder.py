import numpy as np
import math as m
import pandas as pd

def calc_center_of_masses(input_file, mass_dict, snapshots_to_read, num_atoms, atoms_per_adsorbate, MOF_atoms):
    # stores all center of mass values
    COMs = []
    # err counts the number of molecule groupings that do not match the total number of atoms_per_adsorbate.
    err = 0
    #octahedral = [534, 718, 1523, 1538, 1575, 1565, 1643, 1649, 1723, 1777, 1414, 1420, 1430, 1436, 1483, 1493, 1481, 1459, 1405, 1480, 1443, 1432, 1418, 1424, 1505, 1487, 1532, 1549, 1563, 1569, 1586, 1653, 634, 1647, 1627, 1579, 1542, 1701, 1765, 1788, 1727, 1818, 641, 1497, 498, 489, 589, 543, 678, 727, 671, 597, 801, 880, 1398, 1634, 1450, 2383, 2510, 2370, 2387, 2506, 2538, 2714, 2281, 2288, 2293, 2300, 2374, 2429, 2437, 2517, 2522, 2529, 2591, 2596, 3197, 3209, 3216, 3240, 3250, 3258, 3263, 3265, 3298, 3302, 3313, 3336, 3345, 3358, 3380, 3385, 3395, 3417, 3438, 3445, 3469, 3479, 3512, 3543, 3597, 3619, 3633, 2284, 2296, 2305, 2309, 2318, 2363, 2378, 2408, 2418, 2425, 2432, 2453, 2462, 2469, 2490, 2499, 2513, 2525, 2534, 2547, 2587, 2608, 2621, 2692, 3201, 3204, 3213, 3226, 3233, 3236, 3246, 3252, 3262, 3271, 3294, 3309, 3315, 3324, 3334, 3340, 3350, 3354, 3360, 3369, 3378, 3384, 3391, 3397, 3406, 3415, 3421, 3433, 3442, 3455, 3465, 3475, 3516, 3529, 3539, 3593, 3608, 3623, 3638, 3646, 3291, 3426, 3278, 3287, 3430, 3462, 3630]
    octahedralBoyd = [1723, 1777, 1701, 1765, 1788, 1727, 1818, 801, 880, 2714, 2591, 2596, 3512, 3543, 3597, 3619, 3633, 2587, 2608, 2621, 2692, 3516, 3529, 3539, 3593, 3608, 3623, 3638, 3646, 3630]


#    tetrahedral = [2325, 2335, 2341, 2352, 2360, 2383, 2398, 2404, 2443, 2449, 2480, 2486, 2510, 2544, 2554, 2564, 2618, 2628, 2638, 2682, 2688, 2712, 2727, 2735, 2508, 2710, 2716, 2725, 2389, 2370, 2538, 2612, 2714, 2724, 2281, 2288, 2293, 2300, 2315, 2320, 2365, 2374, 2411, 2423, 2429, 2437, 2456, 2467, 2493, 2504, 2517, 2522, 2529, 2549, 2591, 2596, 2603, 2623, 2695, 2706, 2717, 3240, 3250, 3258, 3263, 3265, 3302, 3313, 3358, 3385, 3395, 3469, 3479, 3543, 3553, 3587, 3597, 3625, 3633, 2471, 2515, 2527, 2589, 2601, 2536, 2610, 3640, 2286, 2298, 2307, 2380, 2427, 2434, 2309, 2318, 2363, 2408, 2418, 2453, 2462, 2490, 2499, 2547, 2621, 2692, 2701, 3646, 2719]
    # Iterate through every desired snapshot
    tetrahedralBoyd = [2727, 2735, 2716, 2725, 2714, 2724, 2717, 3625, 3633, 3640, 3646, 2719]

    # Tetra cages are defined by the 3 Zr boned to the mu3-O pointing into the cage. 3 Zr on 4 SBU = 12 total atoms
    tetrahedral = [1815,1821,1803,1819,1808,1820,1812,1822,3633,3625,3638,3630] #mu3OH
  #  tetrahedral = [892, 908, 900, 910, 1815, 1821, 895, 1819, 2726, 2714, 3633, 2717] #mu3O

    #RoggeTetramu3O = [892, 908, 900, 910, 1815, 1821, 895, 1819, 2726, 2714, 3633, 2717]

# Octahedral cage is defined by 5 atoms at each SBU: a central Zr and the 4 bonded linker oxygen. With 6 SBUs forming the corners of the octahedral, this yields 30 atoms
    octahedral = [564, 612, 516, 468, 1500, 1548, 1440, 1819, 708, 1632, 2714, 3276, 3360, 3633, 3420, 660, 2376, 2472, 3637, 2556, 3647, 2400, 2496, 3336, 1608, 2316, 3216, 3630, 3480, 1416]
    for snp in range(snapshots_to_read):
        print('Snapshot Number:', snp + 1)
        # Reads in 1 snapshot at a time
        skips = num_atoms*snp + 9*(snp+1)
        df = pd.read_csv(input_file, skiprows = skips, nrows = num_atoms, error_bad_lines = False, names = ['Atom', 'mol', 'type', 'x', 'y', 'z'], delim_whitespace=True)
        #print(df.head())
        # Create df of tetra and octahedral cell atoms (using inputed lists of atoms for 1 octahedral and 1 tetrahedral cage)
        #df = df[df['mol'].isin(MOF_atoms)]
        dfTetra = df[df['Atom'].isin(tetrahedral)]
        dfOcta = df[df['Atom'].isin(octahedral)]
       # print(len(dfTetra))

        # Map mass values for each element to all atoms in each df
        dfTetra['mass'] = dfTetra['type'].map(MOF_masses)
        dfOcta['mass'] = dfOcta['type'].map(MOF_masses)
        #print(dfTetra.head())
        #print(len(dfTetra), len(dfOcta))

        # Calculate total masses for each cage
        OCTACOMDEN = np.sum(dfOcta['mass'])
        TETRACOMDEN = np.sum(dfTetra['mass'])
       # print(OCTACOMDEN, TETRACOMDEN)

        # Function for calculating the COM in 1 dimension
        # Inputs are atom position list in 1 dimension, masses of each atom, and the total mass of the group of atoms,
        CALCCOMDIR = lambda pos, mass, sumMass : np.sum(pos * mass) / sumMass

        posCOMT = []
        posCOMO = []
        # Calculate the COM in each dimension for each cage
        for pos in ['x', 'y', 'z']:
            posCOMT = posCOMT + [CALCCOMDIR(dfTetra[pos], dfTetra['mass'], TETRACOMDEN)]
            posCOMO = posCOMO + [CALCCOMDIR(dfOcta[pos], dfOcta['mass'], OCTACOMDEN)]
         #   print(len(tetrahedral), len(octahedral))
    return posCOMT, posCOMO


snapshots_to_read = 1 # Read 1 snapshot to calculate COM, could look at more to see how COM changes as framework flexes, but it should remain more or less the same
num_atoms = 4048 # Number of atoms in simulation
atoms_per_adsorbate = 4 #
MOF_atoms = [444] # MOF molecule inteditifier(s)
input_file = r'D:/LAMMPSTraj/RoggeHbond/dump.RoggeTraPPEHbond0Load1.lammpstrj'
#MOF_masses = {'C': 12.01, 'H': 1.01, 'O': 15.999, 'Zr': 91.224} # masses of atoms in the MOF
MOF_masses = {1: 12.01, 2: 12.01, 3: 12.01, 4: 1.01, 5: 1.01, 6: 15.999, 7: 15.999, 8: 15.999, 9: 91.244}
# Call to center of mass calculator
comT, comO = calc_center_of_masses(input_file, MOF_masses, snapshots_to_read, num_atoms, atoms_per_adsorbate, MOF_atoms) # Calculate COM call
print('Cage 1 COM', comT)
print('Cage 2 COM', comO)
difs = [0, 0, 0]
comT = np.array(comT)
comO = np.array(comO)
distCalc = lambda L1,L2 : np.sqrt((L1[0] - L2[0]) ** 2 + (L1[1] - L2[1]) ** 2 + (L1[2] - L2[2]) ** 2)
distval = distCalc(comT, comO)
print('Distance from center to center: ', distval)

difs = comO - comT


print('Distance unit vector: ', difs/distval)



