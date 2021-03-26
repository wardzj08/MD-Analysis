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
UiO66Df = pd.DataFrame({'atom': [1,2,3,4,5],
					    'epsilon': [.10500, .04400, .06000, .06000, .069000],
					    'sigma': [3.430851, 2.571134, 3.118146, 3.118146, 2.783168]
					    })
# Identical df
AceDf = pd.DataFrame({'atom': [ 6,7, 8,9],
                      'epsilon': [0.1569, 0.07948, 0.1947, 0.1947],
                      'sigma': [3.050, 3.820, 3.750, 3.750]
                      })

#IPAdf = pd.DataFrame({'atom': [6, 7, 8, 9, 10],
 #                     'epsilon': [0, 0.1846, .0199, 0.1946, 0.1946],
  #                    'sigma': [0, 3.02, 4.3, 3.75, 3.75]
   #                   })

# Unique identifiers for Atoms in each group
AceAtoms = [6,7,8,9]
UiO66Atoms = [1,2,3,4,5]
IPAAtoms = [6,7,8,9,10]
elementsAce = ['C', 'H', 'O', 'O(mu3)', 'Zr', 'O(Ace)', 'C(Ace)', 'CH3(Ace)', 'CH3(Ace)']
elementsIPA = ['C', 'H', 'O', 'O(mu3)', 'Zr', 'H(IPA)', 'O(IPA)', 'CH(IPA)', 'CH3(IPA)', 'CH3(IPA)']
IPA = False
if IPA == True:
	df1 = pd.concat([UiO66Df, IPAdf], ignore_index = True)
	df2 = df1
	elements = elementsIPA
	AAtoms = IPAAtoms
else:
	df1 = pd.concat([UiO66Df, AceDf], ignore_index = True)
	df2 = df1
	elements = elementsAce
	AAtoms = AceAtoms
#print(len(df1))
#print(len(df2))
# Element for each atom #. This is not necessary, but adds a comment to the end of each
# pair_coeff line, stating the two elements envolved in the parameters on that line


# Create pairs list
#print(df1.index.repeat(9))
DfRep1 = df1[['atom', 'epsilon', 'sigma']].loc[df1.index.repeat(df2.shape[0])].reset_index(drop=True)
DfRep1.columns = ['Atom1', '1eps', '1sigma']
DfRep2 = pd.concat([df2[['atom', 'epsilon', 'sigma']]]*df1.shape[0], ignore_index=True)
DfRep2.columns = ['Atom2', '2eps', '2sigma']
print(len(DfRep1))
print(len(DfRep2))
# Functions to calculate mixed parameters, using arithmatic method
# epsilon_mixed = sqrt(epsilon_i, epsilon_j)
# sigma_mixed = 0.5(sigma_i, sigma_j)
mixedEpsilon = lambda eps1, eps2 : np.sqrt(eps1 * eps2)
mixedSigma = lambda sig1, sig2 : .5 * (sig1 + sig2)

# Preform calculations
mxEp = pd.DataFrame(mixedEpsilon(DfRep2['2eps'], DfRep1['1eps']), columns = ['MixedEpsilon'])
mxSig = pd.DataFrame(mixedSigma(DfRep2['2sigma'], DfRep1['1sigma']), columns = ['MixedSigma'])

# final values
# format: atom1, atom2, mixed epsilon value, mixed sigma value
mixed = pd.concat([DfRep1['Atom1'], DfRep2['Atom2'], mxEp, mxSig], axis=1).reset_index(drop=True)
#print(mixed)

#elements = elementsAce
# Write values to text file
# Write in the pair_coeff format for LAMMPS
write_pairpotentials_tofile = 'DREIDINGmixedljcoul.dat'
with open(write_pairpotentials_tofile,'w') as f:
	for val in mixed.values:
		# If atom#1 is greater than atom#2 it is a repeat of a previous pair, and can be removed
		if val[0] > val[1]:
			continue
		# If atom#1 is different from atom#2, use lj/cut, no charge interactions
		if (val[0] in UiO66Atoms and val[1] in AAtoms) or (val[0] in AAtoms and val[1] in UiO66Atoms):
			f.write(f'pair_coeff {int(val[0])} {int(val[1])} {round(val[2], 6)} {round(val[3], 6)} # {elements[int(val[0])-1]} {elements[int(val[1])-1]}\n')
		# else use lj/cut/coul/long, to include charge interactions
		else:
			f.write(f'pair_coeff {int(val[0])} {int(val[1])} {round(val[2], 6)} {round(val[3], 6)} # {elements[int(val[0])-1]} {elements[int(val[1])-1]}\n')
