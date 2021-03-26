# Converts
import pandas as pd
import numpy as np

COMs = [] # COMs of guest molecules
err = 0 # err counts the number of molecule groupings that do not match the total number of atoms_per_adsorbate.
df = pd.DataFrame(columns = ["Type", "Mol", "x", "y", "z"])
mol_num = 0
mol_prev = 0
atoms_per_mol = 4 # N H H H
masses = {'N': 14.01, 'H': 1.007} # massses associated with each atom

"""Block for reading atoms from a pdb file.
   This will read in all atoms in the file, regardless of being part of seperate models.
   This is good for looking at averages across snapshots
   """
f = './data/NH3_1216.pdb' # pdb file with guest molecules only
for line in open(f):
    v = line.split()
    id = v[0]

    # read in atoms. Need the xyz coords, the atom type (element), and the molecule number
    # Molecules number is not necessarily a field in PDB file. But the atoms should be written in
    # sequential order. So every n atoms are 1 molecules, where n is the number of atoms in the molecule
    if id == 'ATOM':
        print(v)
        type = v[2]
        x = float(v[4])
        y = float(v[5])
        z = float(v[6])
        df = df.append({"Type": type, "Mol": mol_num, "x": x, "y": y, "z":z}, ignore_index=True)
        # if at the final atom of a molecules, increase the molecule identifier by 1
        if int(v[1])%atoms_per_mol == 0:
            mol_num += 1

n = 0 # molecule number tracker
# Sorting of dataframe so entries are grouped by molecule number (each is 1 adsorbate molecule) and sorted
# by element symbol within those groups.
df = df.sort_values('Mol')
df_grouped = df.groupby(df['Mol'])
df_grouped = df_grouped.apply(lambda x: x.sort_values(['Type'], ascending = False))
df_grouped = df_grouped.reset_index(drop = True)
df_grouped = df_grouped.groupby('Mol')

molecule = [] # store molecule positions
for mol, group in df_grouped:
   # at the first molecule, calculate the total mass of the molecule, don't recalculate at additional molecules
    if n == 0:
        # Mass sum, the denominator of the COM equation, calculated during the first molecule run through.
        mass_sum = 0
        for j in range(atoms_per_mol):
                temp_atm = df_grouped.get_group(mol).iloc[j]
                mass_sum += masses[temp_atm['Type']]
        print("Mass of guest molecule: ", mass_sum)
    n += 1 # increase the molecule number

    # Check for groups that have atoms which do not add to the number of atoms in the adsorbate, ideally = 0
    # if error returns a non-zero number, something's wrong
    if len(df_grouped.get_group(mol)) != atoms_per_mol:
            err += 1
            continue

    #Iterate over each atom in the adsorbate, calculating the sum of position * mass over all atoms in the molecule. Add final value to molecule list
    atom_calc = 0
    for k in range(atoms_per_mol):
            atm = df_grouped.get_group(mol).iloc[k]
          #  print('\n',atm)
            atom_calc += masses[atm['Type']] * atm[['x', 'y', 'z']]
    molecule.append(atom_calc)

# Calculate the center of mass using the values in the molecule list and the overall mass
# COM calculation: SUM(Mass of atom n * (xn,yn,zn)) for n = 1 through n = num_adsorbate / mass_sum
COM_SUM = 0
for mol in range(mol_num):
    COM_SUM = molecule[mol]
    COM = COM_SUM / mass_sum
    # Add center of mass to list of COMs
    COMs.append(list(COM))

com_df = pd.DataFrame(np.array(COMs))
com_df.columns = ['x', 'y', 'z']
#com_df.to_csv('COM.csv')
com_df['id'] = list(range(len(com_df['x'])))
com_df['mol'] = [100 for x in range(len(com_df['x']))]
com_df['type'] = [1 for x in range(len(com_df['x']))]

fname_out = 'dump.COMNH3_1216A.lammpstrj'
with open(fname_out, 'w') as f:
    f.write(f'ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{len(COM)}\nITEM: BOX BOUNDS pp pp pp\n0.000e+00 {max(com_df.x)}\n0.000e+00 {max(com_df.y)}\n0.000e+00 {max(com_df.z)}\nITEM: ATOMS ')
com_df.to_csv(fname_out, sep = ' ', columns = ['id', 'mol', 'type', 'x', 'y', 'z'], index = False, mode = 'a')

print('Number of Molecules in Distribution:', len(COMs))
print('Number of Molecules with less than the number of atoms in adsorbate molecule and thus not in the distribution:',err)