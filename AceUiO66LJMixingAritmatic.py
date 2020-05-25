# Compute mixed lj potential parameters
# Uses the arithmatic equations as stated by LAMMPS
# Pass in two dataframes
# Each holding atom type numbers, epsilon values, and sigma values (parameters for lj equation)
# Calculates mixed parameters for every pair of atoms between the two dfs
# If all mixed parameters for a set of atoms is desired, pass identical lists for the two inputs
# This will return all parameters, mixed and self-paired, as the self paired calculations reduce down
# to their starting epsilon and sigma values
# Writes the mixed values to an output file in the lj/coul/cut/long or lj/coul LAMMPS format depending on if the charge interaction should be included:
# pair_coeff ATOM#1 ATOM#2 EPSILONMIX SIGMAMIX [LJcutoff, CoulombicCutoff]
import numpy as np
import pandas as pd

# Set atom numbers and their parameters
# All atoms in the acetone-UiO66 system
UiO66Df = pd.DataFrame({'atom': [1,2,3,4,5, 6, 7, 8, 9],
					    'epsilon': [.10500, .04400, .06000, .06000, .069000, 0.1569, 0.07948, 0.1947, 0.1947],
					    'sigma': [3.430851, 2.571134, 3.118146, 3.118146, 2.783168, 3.050, 3.820, 3.750, 3.750]
					    })
# Identical df
AceDf = pd.DataFrame({'atom': [1,2,3,4,5, 6, 7, 8, 9],
                      'epsilon': [.10500, .04400, .06000, .06000, .069000, 0.1569, 0.07948, 0.1947, 0.1947],
                      'sigma': [3.430851, 2.571134, 3.118146, 3.118146, 2.783168, 3.050, 3.820, 3.750, 3.750]
                      })

# Unique identifiers for Atoms in each group
AceAtoms = [6,7,8,9]
UiO66Atoms = [1,2,3,4,5]

# Element for each atom #. This is not necessary, but adds a comment to the end of each
# pair_coeff line, stating the two elements envolved in the parameters on that line
elements = ['C', 'H', 'O', 'O(mu3)', 'Zr', 'O(Ace)', 'C(Ace)', 'CH3(Ace)', 'CH3(Ace)']

# Create pairs list
UiO66DfRep = UiO66Df[['atom', 'epsilon', 'sigma']].loc[UiO66Df.index.repeat(AceDf.shape[0])].reset_index(drop=True)
UiO66DfRep.columns = ['Atom1', '1eps', '1sigma']
AceDfRep = pd.concat([AceDf[['atom', 'epsilon', 'sigma']]]*UiO66Df.shape[0], ignore_index=True)
AceDfRep.columns = ['Atom2', '2eps', '2sigma']

# Functions to calculate mixed parameters, using arithmatic method
# epsilon_mixed = sqrt(epsilon_i, epsilon_j)
# sigma_mixed = 0.5(sigma_i, sigma_j)
mixedEpsilon = lambda eps1, eps2 : np.sqrt(eps1 * eps2)
mixedSigma = lambda sig1, sig2 : .5 * (sig1 + sig2)

# Preform calculations
mxEp = pd.DataFrame(mixedEpsilon(AceDfRep['2eps'], UiO66DfRep['1eps']), columns = ['MixedEpsilon'])
mxSig = pd.DataFrame(mixedSigma(AceDfRep['2sigma'], UiO66DfRep['1sigma']), columns = ['MixedSigma'])

# final values
# format: atom1, atom2, mixed epsilon value, mixed sigma value
mixed = pd.concat([UiO66DfRep['Atom1'], AceDfRep['Atom2'], mxEp, mxSig], axis=1).reset_index(drop=True)
print(mixed)

# Write values to text file
# Write in the pair_coeff format for LAMMPS
write_pairpotentials_tofile = 'mixedLJCOULPotentials.dat'
with open(write_pairpotentials_tofile,'w') as f:
	for val in mixed.values:
		# If atom#1 is greater than atom#2 it is a repeat of a previous pair, and can be removed
		if val[0] > val[1]:
			continue
		# If atom#1 is different from atom#2, use lj/cut, no charge interactions
		if (val[0] in UiO66Atoms and val[1] in AceAtoms) or (val[0] in AceAtoms and val[1] in UiO66Atoms):
			f.write(f'pair_coeff {int(val[0])} {int(val[1])} lj/cut {round(val[2], 6)} {round(val[3], 6)} # {elements[int(val[0])-1]} {elements[int(val[1])-1]}\n')
		# else use lj/cut/coul/long, to include charge interactions
		else:
			f.write(f'pair_coeff {int(val[0])} {int(val[1])} lj/cut/coul/long {round(val[2], 6)} {round(val[3], 6)} # {elements[int(val[0])-1]} {elements[int(val[1])-1]}\n')
