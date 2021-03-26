import pandas as pd
import numpy as np

dfAce = pd.DataFrame([[5.3, 15.8, 13.571], [5.3, 15.8, 14.8,], [5.3, 17.097397, 15.591934], [5.3, 14.502603, 15.591934], columns = ['x', 'y', 'z'])


dfAce['mass'] = [15.999, 12.009, 15.0345, 15.0345]
print(dfAce)

AceCOMDEN = np.sum(dfAce['mass'])
  #      print(OCTACOMDEN, TETRACOMDEN)

        # Function for calculating the COM in 1 dimension
        # Inputs are atom position list in 1 dimension, masses of each atom, and the total mass of the group of atoms,
CALCCOMDIR = lambda pos, mass, sumMass : np.sum(pos * mass) / sumMass

        # Calculate the COM in each dimension for each cage
for pos in ['x', 'y', 'z']:
    posCOMT = CALCCOMDIR(dfAce[pos], dfAce['mass'], AceCOMDEN)
    print(posCOMT)