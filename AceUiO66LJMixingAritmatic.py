# Compute mixed lj potential parameters
# Uses the aritmatic equations as stated in LAMMPS
import numpy as np
import pandas as pd

# Set atom numbers and their parameters
UiO66Df = pd.DataFrame({'atom': [1,2,3,4,5], 'epsilon': [.10500, .04400, .06000, .06000, .069000], 'sigma': [3.430851, 2.571134, 3.118146, 3.118146, 2.783168]})
AceDf = pd.DataFrame({'atom': [6, 7, 8, 9], 'epsilon': [0.1569, 0.07948, 0.1947, 0.1947], 'sigma': [3.050, 3.820, 3.750, 3.750]})


# Create a pairs list
UiO66DfRep = pd.concat([UiO66Df[['atom', 'epsilon', 'sigma']]]*AceDf.shape[0], ignore_index=True)
UiO66DfRep.columns = ['Atom1', '1eps', '1sigma']
AceDfRep = AceDf[['atom', 'epsilon', 'sigma']].loc[AceDf.index.repeat(UiO66Df.shape[0])].reset_index(drop=True)
AceDfRep.columns = ['Atom2', '2eps', '2sigma']

# Functions to calculate mixed parameters
# epsilon_mixed = sqrt(epsilon_i, epsilon_j)
# sigma_mixed = 0.5(sigma_i, sigma_j)
mixedEpsilon = lambda eps1, eps2 : np.sqrt(eps1 * eps2)
mixedSigma = lambda sig1, sig2 : .5 * (sig1 +sig2)

# Preform calculations
mxEp = pd.DataFrame(mixedEpsilon(AceDfRep['2eps'], UiO66DfRep['1eps']), columns = ['MixedEpsilon'])
mxSig = pd.DataFrame(mixedSigma(AceDfRep['2sigma'], UiO66DfRep['1sigma']), columns = ['MixedSigma'])

# final values
mixed = pd.concat([UiO66DfRep['Atom1'], AceDfRep['Atom2'], mxEp, mxSig], axis=1).reset_index(drop=True)
print(mixed)

# Write values to text file
# Write in the pair_coeff format for LAMMPS
write_pairpotentials_tofile = 'mixedLJCOULPotentials.dat'
with open(write_pairpotentials_tofile,'w') as f:
	for val in mixed.values:
		f.write(f'pair_coeff {int(val[0])} {int(val[1])} {val[2]} {val[3]} 12.5 0.0\n')
