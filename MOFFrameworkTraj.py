# Given a lammpstrj file of a adsorbate diffusing through a MOF (or really any MOF system) returns a lammpstrj file of the MOF by itself, at timestep = 0; Useful in conjunction with center of mass distribution calculations of the adsorbate, for vizualizing where in the MOF the adsorbate is on average. 

import pandas as pd
import numpy as np
import math as m

input_file = 'dump.production.lammpstrj'
num_atoms = 4048
MOF_mol_num = 444

# Input:
# input_file - name of file to read from
# num_atoms - number of atoms in a single timetstep/snapshot
# MOF_mol_num - molecule number corresponding to the MOF
# output_file - file to write to 
def create_framework_file(input_file, num_atoms, MOF_mol_num, output_file = 'dump.rigidFramework.lammpstrj'):
    df = pd.read_csv(input_file, skiprows = 9, nrows = num_atoms, error_bad_lines = False, names = ['id', 'mol','type', 'x', 'y', 'z'], delim_whitespace=True)
        
    # Removes all atoms that belong to the MOF (in UiO-66 potential framework, mol = 444)
    df = df[df['mol'] == MOF_mol_num]
#df = pd.DataFrame(np.array(coms))
#df.columns = ['x', 'y', 'z']
#df.to_csv('COM.csv')
#df['id'] = list(range(5000, 5000 + len(df['x'])))
#df['mol'] = [444 for x in range(len(df['x']))]
#df['type'] = [2 for x in range(len(df['x']))]

    with open(output_file, 'w') as f:
        f.write('ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n3648\nITEM: BOX BOUNDS pp pp pp\n0.000e+00 4.14008e+01\n0.000e+00 4.14008e+01\n0.000e+00 4.14008e+01\nITEM: ATOMS ')

    df.to_csv(output_file, sep = ' ', index = False, mode = 'a')

create_framework_file(input_file, num_atoms, MOF_mol_num)