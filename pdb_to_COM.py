
import pandas as pd
import numpy as np
COMs = []
    # err counts the number of molecule groupings that do not match the total number of atoms_per_adsorbate.
err = 0
df = pd.DataFrame(columns = ["Type", "Mol", "x", "y", "z"])
mol_num = 0
mol_prev = 0
atoms_per_mol = 4
f = './data/NH3_1216.pdb'
for line in open(f):
    v = line.split()
    id = v[0]

    # Only need to read in atoms. Need the xyz coords, the atom type (element), and the molecules number
    if id == 'ATOM':
        print(v)
        type = v[2]
        x = float(v[4])
        y = float(v[5])
        z = float(v[6])
        df = df.append({"Type": type, "Mol": mol_num, "x": x, "y": y, "z":z}, ignore_index=True)
        # if at the final atom of a molecules, increase the molecule identifier
        if int(v[1])%atoms_per_mol == 0:
            print("True")
            mol_num += 1

            #if(mol_num) < mol_prev:
             #   break
            #mol_prev = mol_num
   #    if mol_num >= 1000:
    #        break
#f.close()
#print(df)
   #     type = v[2]
    #    if type == 'CA':
     #       residue = v[3]
      #      type_of_chain = v[4]
       #     atom_count = int(v[5])
        #    position = v[6:8]

masses = {'N': 14.01, 'H': 1.007}
n = 0
        # Removes all atoms that belong to the MOF (in UiO-66 potential framework, mol = 444)
        # If Framework has more than one mol number, passing a list of mol values will work
df = df.sort_values('Mol')
#print(df)
        # Sorting of dataframe so entries are grouped by molecule number (each is 1 adsorbate molecule) and sorted
        # by element symbol within those groups.
df_grouped = df.groupby(df['Mol'])
#print(df_grouped)
df_grouped = df_grouped.apply(lambda x: x.sort_values(['Type'], ascending = False))
df_grouped = df_grouped.reset_index(drop = True)
df_grouped = df_grouped.groupby('Mol')

#print(df_grouped)
molecule = []
for mol, group in df_grouped:
   # print(group)
    if n == 0:

                # Mass sum, the denominator of the COM equation, calculated during the first molecule run through.
                # Becuase the atoms are grouped by mol number and sorted by element, this is identical for all molecules
        # and only needs to be computed once
        mass_sum = 0
        for j in range(atoms_per_mol):
                    temp_atm = df_grouped.get_group(mol).iloc[j]
                    #print(temp_atm)
                    mass_sum += masses[temp_atm['Type']]
                    print(mass_sum)

    n += 1
            # Check for groups that have atoms which do not add to the number of atoms in the adsorbate, ideally = 0
    if len(df_grouped.get_group(mol)) != atoms_per_mol:
            err += 1
            continue



            #Iterate over each atom in the adsorbate, calculating the position * mass. Add to the molecule list
    atom_calc = 0
    for k in range(atoms_per_mol):
            atm = df_grouped.get_group(mol).iloc[k]
          #  print('\n',atm)
            atom_calc += masses[atm['Type']] * atm[['x', 'y', 'z']]
    molecule.append(atom_calc)

            # Calculate the center of mass using the values in the molecule list and the overall mass
            # COM calculation: SUM(Mass of atom n * (xn,yn,zn)) for n = 1 through n = num_adsorbate / mass_sum

COM_SUM = 0
#print(molecule)
for mol in range(mol_num):
    #for atom in range(atoms_per_mol):
    COM_SUM = molecule[mol]
    COM = COM_SUM / mass_sum
            #print(mol)
            #print(COM)
            # Add center of mass to list of COMs
    COMs.append(list(COM))


print(COMs)
print(len(COMs))


com_df = pd.DataFrame(np.array(COMs))
com_df.columns = ['x', 'y', 'z']
#com_df.to_csv('COM.csv')
com_df['id'] = list(range(len(com_df['x'])))
com_df['mol'] = [100 for x in range(len(com_df['x']))]
com_df['type'] = [1 for x in range(len(com_df['x']))]

with open('dump.COMNH3_1216.lammpstrj', 'w') as f:
    f.write('ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n80000\nITEM: BOX BOUNDS pp pp pp\n0.000e+00 70.14008e+01\n0.000e+00 70.14008e+01\n0.000e+00 70.14008e+01\nITEM: ATOMS ')
com_df.to_csv('dump.COMNH3_1216.lammpstrj', sep = ' ', columns = ['id', 'mol', 'type', 'x', 'y', 'z'], index = False, mode = 'a')

print('Number of Molecules in Distribution:', len(COMs))
print('Number of Molecules with less than the number of atoms in adsorbate molecule and thus not in the distribution:',err)