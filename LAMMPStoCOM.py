# Converts
import pandas as pd
import numpy as np

COMs = [] # COMs of guest molecules
err = 0 # err counts the number of molecule groupings that do not match the total number of atoms_per_adsorbate.
df = pd.DataFrame(columns = ["Type", "Mol", "x", "y", "z"])
loading = 3.125/2
mol_num = int(loading*32)
mol_prev = 0
atoms_per_mol = 5 # N H H H
masses = {10: 1.0079, 11: 15.9994, 12: 13.0186, 13: 15.0345, 14: 15.0345} # massses associated with each atom
snapshots_to_read = 1
MOF_atoms = [1]
input_file = 'data/dump.0LoadIPA.lammpstrj'
num_atoms = 3648 + mol_num*5

"""Block for reading atoms from a pdb file.
   This will read in all atoms in the file, regardless of being part of seperate models.
   This is good for looking at averages across snapshots
   """
for snp in range(snapshots_to_read):
        print('Snapshot Number:', snp + 1)
        # Reads in 1 snapshot at a time
        if snp == 0:
           df = pd.read_csv(input_file, skiprows = 9, nrows = num_atoms, names = ['mol', 'type', 'x', 'y', 'z'], usecols=[1,2,3,4,5], delim_whitespace=True)
        else:
            skips = num_atoms*snp + 9*(snp+1)
            df = pd.read_csv(input_file, skiprows = skips, nrows = num_atoms, error_bad_lines = False, names = ['mol', 'type', 'x', 'y', 'z'], usecols=[1,2,3,4,5], delim_whitespace=True)
        #print(df.head())
        # Removes all atoms that belong to the MOF (in UiO-66 potential framework, mol = 444)
        # If Framework has more than one mol number, passing a list of mol values will work
        df = df[~df['mol'].isin(MOF_atoms)].sort_values('mol')
        print(df.head(10))
        print(len(df))
        n = 0 # molecule number tracker
        # Sorting of dataframe so entries are grouped by molecule number (each is 1 adsorbate molecule) and sorted
        # by element symbol within those groups.
        df = df.sort_values('mol')
        df_grouped = df.groupby(df['mol'])
        df_grouped = df_grouped.apply(lambda x: x.sort_values(['type'], ascending = False))
        df_grouped = df_grouped.reset_index(drop = True)
        df_grouped = df_grouped.groupby('mol')

        molecule = [] # store molecule positions
        for mol, group in df_grouped:
   # at the first molecule, calculate the total mass of the molecule, don't recalculate at additional molecules
            if n == 0:
        # Mass sum, the denominator of the COM equation, calculated during the first molecule run through.
                mass_sum = 0
                for j in range(atoms_per_mol):
                    temp_atm = df_grouped.get_group(mol).iloc[j]
                    mass_sum += masses[temp_atm['type']]
            #print("Mass of guest molecule: ", mass_sum)
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
                atom_calc += masses[atm['type']] * atm[['x', 'y', 'z']]
            molecule.append(atom_calc)
            #print("Atom calc", atom_calc)
            #print("Mol", atm)
            print("Mol", molecule)
# Calculate the center of mass using the values in the molecule list and the overall mass
# COM calculation: SUM(Mass of atom n * (xn,yn,zn)) for n = 1 through n = num_adsorbate / mass_sum

        COM_SUM = 0
        for mol in range(mol_num):
            COM_SUM = molecule[mol]
            COM = COM_SUM / mass_sum

            #print("TEST",COM_SUM)
    # Add center of mass to list of COMs
            COMs.append(list(COM))
        print(COMs)
com_df = pd.DataFrame(np.array(COMs))
print(com_df.head())
print(len(com_df))
com_df.columns = ['x', 'y', 'z']
#com_df.to_csv('COM.csv')
com_df['id'] = list(range(len(com_df['x'])))
com_df['mol'] = [100 for x in range(len(com_df['x']))]
com_df['type'] = [1 for x in range(len(com_df['x']))]

fname_out = 'dump.k.lammpstrj'
with open(fname_out, 'w') as f:
    f.write(f'ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{len(com_df)}\nITEM: BOX BOUNDS pp pp pp\n0.000e+00 {max(com_df.x)}\n0.000e+00 {max(com_df.y)}\n0.000e+00 {max(com_df.z)}\nITEM: ATOMS ')
com_df.to_csv(fname_out, sep = ' ', columns = ['id', 'mol', 'type', 'x', 'y', 'z'], index = False, mode = 'a')

print('Number of Molecules in Distribution:', len(COMs))
print('Number of Molecules with less than the number of atoms in adsorbate molecule and thus not in the distribution:',err)