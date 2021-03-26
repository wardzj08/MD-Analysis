import pandas as pd
import numpy as np
import math
# Shows all columns of dataframe instead of the front and end only
pd.set_option('display.max_columns', None)

input_file = 'dump.IPAUnwrapped' # file to read trajectory from
snapshots_to_read = 500# Number of timesteps to read from traj file

loading = 3.125 # level of loading  (diffusing molecules per primative cell), 3.125 = 0 loading
atoms_per_adsorbate = 5 # number of atoms in the diffusing molecule 5 atoms in IPA


num_atoms = int(32 * loading * atoms_per_adsorbate)
print('Atoms in system: ', num_atoms)

num_mols = int(num_atoms / atoms_per_adsorbate)
print('Molecules in system:', num_mols)
mass_dict = {10: 1.0079, 11: 15.9994, 12: 13.0186, 13: 15.0345, 14: 15.0345}
MolMass = sum(mass_dict.values())
print("Mass of adsorbate molecule:", MolMass)

#Define box size
# Since coordinates are shifted to a 0 minimum (so 0,0,0 is a corner of the box), only need the length of the box
BoxSize = 41.9567985

# Calculate x1 - x0, where x1 is the position (x,y or z) at the current time, and x0 is initial position.
def distance(x1, x0):
    x0 = x0 - np.floor(x0/BoxSize) * BoxSize
    dx = x1 - x0
    dx = dx - np.round(dx/BoxSize) * BoxSize
    return dx
fname1 = "TSTMSD3np"
with open(fname1, "w") as out:
	out.write("1\n2\n3\n")
start = 4000
# loop over each snapshot
#msd_ds = np.zeros(shape = (snapshots_to_read, 4))
#msd_dsat = np.zeros(shape = (snapshots_to_read, 4))
#msd_dc = np.zeros(shape = (snapshots_to_read, 4))
for snp in [0] + list(range(start, snapshots_to_read)):
	if snp%10 == 0:

		print('\nSnapshot Number:', snp + 1)
        # Reads in 1 snapshot at a time
        skips = num_atoms*snp + 9*(snp+1)
        df = pd.read_csv(input_file, skiprows = skips, nrows = num_atoms, error_bad_lines = False, names = ['Atom','mol', 'type', 'x', 'y', 'z'], usecols = ['Atom', 'type', 'x', 'y', 'z'], delim_whitespace=True)

	# sort by mol number, keeps the order consistant from timestep to timstep
        df = df.sort_values('Atom')
        # Add masses to df
        df['Mass'] = df['type'].map(mass_dict)
        df = df.values
	#print(df.head(6))
        # shift to a 0 minimum
        #print(df.head())
        #df[['x', 'y', 'z']] += 10.83147275
        #print(df.head())
        # calculate COM of each molecule (x,y,z)
        # store COM data here
        #df_mol = pd.DataFrame(columns = ['Mol', 'x', 'y', 'z'])
        # Calculate the COM in each dimension for each molecule
	"""        for mol in df['mol'].unique():
            COM_temp = {'Mol': int(mol), 'x': 0, 'y': 0, 'z': 0}
            df_1mol = df[df['mol'] == mol].sort_values('type')
            #print(df_1mol)
            # Find COM for each direction
            for pos in ['x', 'y', 'z']:
                pos1 = df_1mol[pos]
                # shift to a 0 minimum
                #pos1 +=  10.83147275
                pos1 = pos1/BoxSize * 2 * math.pi
                cosPos1 = np.cos(pos1)
                sinPos1 = np.sin(pos1)
                #print(pos1)
                #print(cosPos1)
                #print(df_1mol)
                cosAvg = np.sum(cosPos1 * df_1mol['Mass']) / MolMass
                sinAvg = np.sum(sinPos1 * df_1mol['Mass']) / MolMass
                #print(cosAvg)
                atan = math.atan2(-sinAvg, -cosAvg) + math.pi
                posCOM = BoxSize * atan/(2*math.pi)
              #  print(posCOM)
                COM_temp[pos] = posCOM
           # print(COM_tmp)
           # add to the molecule COM list
            df_mol = df_mol.append(COM_temp, ignore_index=True)
	    """
	"""
        COM_sys = {'x': 0, 'y': 0, 'z': 0}
        df_1 = df
            #print(df_1mol)
            # Find COM for each direction
        for pos in ['x', 'y', 'z']:
            pos1 = df_1[pos]
                # shift to a 0 minimum
                #pos1 +=  10.83147275
            pos1 = pos1/BoxSize * 2 * math.pi
            cosPos1 = np.cos(pos1)
            sinPos1 = np.sin(pos1)
                #print(pos1)
                #print(cosPos1)
                #print(df_1mol)
            cosAvg = np.sum(cosPos1 * df_1['Mass']) / MolMass*num_mols
            sinAvg = np.sum(sinPos1 * df_1['Mass']) / MolMass*num_mols
                #print(cosAvg)
            atan = math.atan2(-sinAvg, -cosAvg) + math.pi
            posCOM = BoxSize * atan/(2*math.pi)
              #  print(posCOM)
            COM_sys[pos] = posCOM
           # print(COM_tmp)
           # add to the molecule COM list
       # df_mol = df_mol.append(COM_temp, ignore_index=True)
        # Store initial Molecule COMs
	"""
        if snp == 0:
            df_atom_init = df
        #    df_mol_init = df_mol
      #      COM_init = COM_sys
        # Calculate x-x0, y-y0, z-z0 for each mol COM
	""" diffs = np.zeros(shape = (num_mols, 3))
        for mol in range(num_mols):

            # check if they match
            #print(df_mol.iloc[[mol]])
            #print(df_mol_init.iloc[[mol]])
            xdiff = distance(df_mol.iloc[[mol]]['x'].values[0], df_mol_init.iloc[[mol]]['x'].values[0])# df_mol.iloc[[mol]]['x'] - df_mol_init.iloc[[mol]]['x']
            ydiff = distance(df_mol.iloc[[mol]]['y'].values[0], df_mol_init.iloc[[mol]]['y'].values[0]) #df_mol.iloc[[mol]]['y'] - df_mol_init.iloc[[mol]]['y']
            zdiff = distance(df_mol.iloc[[mol]]['z'].values[0], df_mol_init.iloc[[mol]]['z'].values[0])#df_mol.iloc[[mol]]['z'] - df_mol_init.iloc[[mol]]['z']
           # xdiff = distance(df_mol.iloc[[mol]]['x'].values[0] - COM_sys['x'], df_mol_init.iloc[[mol]]['x'].values[0] - COM_init['x'])# df_mol.iloc[[mol]]['x'] - df_mol_init.iloc[[mol]]['x']
           # ydiff = distance(df_mol.iloc[[mol]]['y'].values[0]- COM_sys['y'], df_mol_init.iloc[[mol]]['y'].values[0] - COM_init['y']) #df_mol.iloc[[mol]]['y'] - df_mol_init.iloc[[mol]]['y']
           # zdiff = distance(df_mol.iloc[[mol]]['z'].values[0] - COM_sys['z'], df_mol_init.iloc[[mol]]['z'].values[0] - COM_init['z'])#df_mol.iloc[[mol]]['z'] - df_mol_init.iloc[[mol]]['z']

            diffs[mol,0] = xdiff
            diffs[mol,1] = ydiff
            diffs[mol,2] = zdiff
        #print("DIFFS")
        #print(diffs)
	"""
	"""     # 1. square each COM and sum in each direction
        #print("SQUARES")
        squares = np.square(diffs)
        #print(squares[1:5, :])
        #print("SUM")
        sum_direction = squares.sum(axis = 0)
       # print(sum_direction)
        MSD = sum_direction / num_mols
        # add total msd
        MSD = np.append(MSD, sum(MSD))
        #print("MSD for Ds at time:", snp*100)
       # print(num_mols, MSD)
        msd_ds[snp, :] = MSD
	"""

 	# Calculate x-x0, y-y0, z-z0 for each mol COM
        diffsat = np.zeros(shape = (num_atoms, 3))
        #for mol in range(num_atoms):

            # check if they match
            #print(df_mol.iloc[[mol]])
            #print(df_mol_init.iloc[[mol]])
	xdiff = distance(df[:, 2], df_atom_init[:, 2])# df_mol.iloc[[mol]]['x'] - df_mol_init.iloc[[mol]]['x']
	ydiff = distance(df[:, 3], df_atom_init[:, 3]) #df_mol.iloc[[mol]]['y'] - df_mol_init.iloc[[mol]]['y']
	zdiff = distance(df[:, 4], df_atom_init[:, 4])#df_mol.iloc[[mol]]['z'] - df_mol_init.iloc[[mol]]['z']
           # xdiff = distance(df.iloc[[mol]]['x'].values[0] - COM_sys['x'], df_atom_init.iloc[[mol]]['x'].values[0]-COM_init['x'])# df_mol.iloc[[mol]]['x'] - df_mol_init.iloc[[mol]]['x']
            #ydiff = distance(df.iloc[[mol]]['y'].values[0] - COM_sys['y'], df_atom_init.iloc[[mol]]['y'].values[0]-COM_init['y']) #df_mol.iloc[[mol]]['y'] - df_mol_init.iloc[[mol]]['y']
            #zdiff = distance(df.iloc[[mol]]['z'].values[0] - COM_sys['z'], df_atom_init.iloc[[mol]]['z'].values[0]-COM_init['z'])#df_mol.iloc[[mol]]['z'] - df_mol_init.iloc[[mol]]['z']

	diffsat[:,0] = xdiff
	diffsat[:,1] = ydiff
        diffsat[:,2] = zdiff
        #print("DIFFS")
        #print(diffs)


        # 1. square each COM and sum in each direction
        #print("SQUARES")
        squaresat = np.square(diffsat)
        #print(squares[1:5, :])
        #print("SUM")
        sum_directionat = squaresat.sum(axis = 0)
       # print(sum_direction)
        MSDat = sum_directionat / num_atoms
        # add total msd
        MSDat = np.append(MSDat, sum(MSDat))
        #print("MSD for Ds at time:", snp*100)
       # print(num_mols, MSD)
    #    msd_dsat[snp, :] = MSDat
	with open(fname1, "a") as out:
		out.write("{} 4".format(100*(snp)))
		out.write("\n")
		out.write("1 {}".format(MSDat[0]))
		out.write("\n")
		out.write("2 {}".format(MSDat[1]))
		out.write("\n")
		out.write("3 {}".format(MSDat[2]))
		out.write("\n")
		out.write("4 {}".format(MSDat[3]))
		out.write("\n")


	"""     # 2. Sum each COM and square in each direction
        sum_diffs = diffs.sum(axis = 0)
        square_sums = np.square(sum_diffs)
        #print("MSD Equivelent for Corrected Diffusivity")
        MSD2 = square_sums / num_mols
        MSD2 = np.append(MSD2, sum(MSD2))
        msd_dc[snp, :] = MSD2
        #print(MSD2)
	"""


#np.savetxt('MSD1TST3.csv', msd_ds, delimiter = ",")
#np.savetxt('MSD12TST3.csv', msd_dsat, delimiter = ",")
#np.savetxt('MSD2TST.csv', msd_dc, delimiter = ",")
