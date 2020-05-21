import numpy as np
import math as m
import pandas as pd

def calc_center_of_masses(input_file, mass_dict, snapshots_to_read, num_atoms, atoms_per_adsorbate, MOF_atoms):
    # stores all center of mass values
    COMs = []
    # err counts the number of molecule groupings that do not match the total number of atoms_per_adsorbate.
    err = 0
    octahedral = [755, 824, 844, 792, 828, 836, 909, 534, 718, 756, 303, 387, 1279, 1523, 1538, 1575, 1565, 1643, 1658, 1649, 1672, 1717, 1723, 1747, 1761, 1777, 1767, 1820, 1816, 821, 827, 1414, 1420, 1430, 1436, 1483, 1493, 1481, 1459, 1405, 1480, 1443, 1432, 1418, 1424, 1505, 1487, 1532, 1549, 1563, 1569, 1586, 1653, 634, 1647, 1627, 1579, 1542, 1701, 1735, 1765, 1788, 1727, 1745, 1811, 1781, 1818, 1721, 641, 1497, 1803, 498, 489, 589, 543, 1771, 678, 727, 671, 759, 748, 597, 801, 841, 843, 873, 880, 895, 1398, 1667, 1736, 1756, 1634, 1668, 1708, 1740, 1748, 1821, 1450, 2383, 2510, 2578, 2644, 2672, 2712, 2727, 2370, 2387, 2506, 2538, 2577, 2612, 2656, 2664, 2708, 2714, 2722, 2724, 2187, 3038, 3127, 3182, 2281, 2288, 2293, 2300, 2374, 2429, 2437, 2517, 2522, 2529, 2569, 2581, 2591, 2596, 2603, 2663, 2669, 2717, 3197, 3209, 3216, 3240, 3250, 3258, 3263, 3265, 3298, 3302, 3313, 3336, 3345, 3358, 3380, 3385, 3395, 3417, 3438, 3445, 3469, 3479, 3512, 3519, 3543, 3553, 3561, 3567, 3597, 3603, 3619, 3625, 3633, 3648, 2645, 3482, 3585, 3640, 2284, 2296, 2305, 2309, 2318, 2363, 2378, 2408, 2418, 2425, 2432, 2453, 2462, 2469, 2490, 2499, 2513, 2525, 2534, 2547, 2587, 2599, 2608, 2621, 2643, 2653, 2692, 2701, 3201, 3204, 3213, 3226, 3233, 3236, 3246, 3252, 3262, 3271, 3294, 3309, 3315, 3324, 3334, 3340, 3350, 3354, 3360, 3369, 3378, 3384, 3391, 3397, 3406, 3415, 3421, 3433, 3442, 3455, 3465, 3475, 3483, 3494, 3507, 3516, 3529, 3539, 3549, 3573, 3583, 3593, 3599, 3608, 3617, 3623, 3638, 3646, 3569, 3642, 3627, 3635, 2572, 2665, 2719, 3291, 3426, 3490, 3500, 3540, 3556, 3584, 3639, 3278, 3287, 3430, 3462, 3489, 3504, 3536, 3568, 3576, 3630, 3634, 3636]

    tetrahedral = [2325, 2335, 2341, 2352, 2360, 2383, 2398, 2404, 2443, 2449, 2480, 2486, 2510, 2544, 2554, 2564, 2618, 2628, 2638, 2682, 2688, 2712, 2727, 2735, 2508, 2710, 2716, 2725, 2389, 2370, 2538, 2612, 2714, 2724, 2281, 2288, 2293, 2300, 2315, 2320, 2365, 2374, 2411, 2423, 2429, 2437, 2456, 2467, 2493, 2504, 2517, 2522, 2529, 2549, 2591, 2596, 2603, 2623, 2695, 2706, 2717, 3240, 3250, 3258, 3263, 3265, 3302, 3313, 3358, 3385, 3395, 3469, 3479, 3543, 3553, 3587, 3597, 3625, 3633, 2471, 2515, 2527, 2589, 2601, 2536, 2610, 3640, 2286, 2298, 2307, 2380, 2427, 2434, 2309, 2318, 2363, 2408, 2418, 2453, 2462, 2490, 2499, 2547, 2621, 2692, 2701, 3646, 2719]
    # Iterate through every desired snapshot
    for snp in range(snapshots_to_read):
        print('Snapshot Number:', snp + 1)
        # Reads in 1 snapshot at a time
        skips = num_atoms*snp + 9*(snp+1)
        df = pd.read_csv(input_file, skiprows = skips, nrows = num_atoms, error_bad_lines = False, names = ['Atom','mol', 'type', 'element', 'x', 'y', 'z'], delim_whitespace=True)

        # Removes all atoms that belong to the MOF (in UiO-66 potential framework, mol = 444)
        # If Framework has more than one mol number, passing a list of mol values will work
        df = df[df['mol'].isin(MOF_atoms)]
        dfTetra = df[df['Atom'].isin(tetrahedral)]
        dfOcta = df[df['Atom'].isin(octahedral)]
#        print(len(dfTetra))

        dfTetra['mass'] = dfTetra['element'].map(MOF_masses)
        dfOcta['mass'] = dfOcta['element'].map(MOF_masses)
 #       print(dfTetra.head())

        OCTACOMDEN = np.sum(dfOcta['mass'])
        TETRACOMDEN = np.sum(dfTetra['mass'])
  #      print(OCTACOMDEN, TETRACOMDEN)

        CALCCOMDIR = lambda pos,mass, sumMass : np.sum(pos * mass) / sumMass

        for pos in ['x', 'y', 'z']:
            posCOMT = CALCCOMDIR(dfTetra[pos], dfTetra['mass'], TETRACOMDEN)
            posCOMO = CALCCOMDIR(dfOcta[pos], dfOcta['mass'], OCTACOMDEN)
   #         print(posCOMT, posCOMO)
    return COMs

# Mass dictionary parameter. Store each unique mass value with its corresponding element, to be called in the COM fucnction
# These are parameters to change based on each system
#masses = {'C': 12.01, 'O': 15.999}#, 'H': 1.007}
snapshots_to_read = 1
num_atoms = 4288
atoms_per_adsorbate = 4
MOF_atoms = [444]
input_file = './data/dump.productionAce5load.lammpstrj'
MOF_masses = {'C': 12.01, 'H': 1.01, 'O': 15.999, 'Zr': 91.224}
# Call to center of mass calculator
coms = calc_center_of_masses(input_file, MOF_masses, snapshots_to_read, num_atoms, atoms_per_adsorbate, MOF_atoms)