# Find position of acetone molecule relative to 3 zirconium and 1 mu3H for hydrogen bonding to engage
# Distances for acetone atoms to the zr,H atoms given from Jon's calculations which show H bonding occurs
# between acetone oxygen and mu3H at 1.8 Angstroms
import numpy as np
import itertools
from scipy import optimize
# MOF ORDER: H H H H Zr1 Zr1 Zr1 Zr1 Zr2 Zr2 Zr2 Zr2 Zr3 Zr3 Zr3 Zr3
# Ace Order: O C CH31 CH32 *4
# Ace-Ace distances added to end in O-C O-CH31 O-CH32 C-O, C-CH31 C-CH32 CH31-O, CH31-Cx, CH31-CH32x, CH32-O, CH32-C, CH32-CH31
dist = [1.8, 2.94, 4.17, 3.64, 4.27, 5.43, 6.48, 6.14, 4.25, 5.36, 6.29, 6.16, 4.37, 5.2, 6.642, 5.17]
AceInterDist = [1.229, 2.40154392, 2.40154392, 1.229, 1.51999949, 1.51999949, 2.40154392, 1.51999949, 2.594794, 2.40154392, 1.51999949, 2.594794]
dist = dist + AceInterDist
#dist = np.zeros(3)
# O C CH31 CH32
xBead = [0, 0, 0, 0]
yBead = [0, 0, 1.297, -1.297]
zBead = [0, 1.229, 2.02, 2.02]

# Pull these from d
# H(atom # 2581) Zr1 (atom # 2717) Zr2 (atom # 3633) Zr3 (atom # 3625)
xMOF = [21.96126, 20.70040, 20.70040, 23.18217]
yMOF = [18.75353, 18.21863, 20.70040, 20.70040]
zMOF = [22.63941, 20.70040, 23.18217, 20.70040]
tstdist = [1, 1]
LocAce = np.transpose(np.array((xBead, yBead, zBead))) # Ace
LocMOF = np.transpose(np.array((xMOF, yMOF, zMOF))) # MOF

# Function to minimize (difference of two distances)
def fct(LOCAce):
	# reshape to xyz coordinates for each atom instead of a list of all coordinates
   LOCAce = np.reshape(LOCAce, (4,3)) # Ace

   # array of all acetone and MOF reference atoms
   LOCS = np.concatenate((LocMOF, LOCAce))
   # calculate cartesian product of Acetone atoms and MOF atoms + Acetone atoms
   CartProd = np.array(list(itertools.product(*[LOCS, LOCAce])))
   # Remove coordinate sets that are two positions for the same atom
   cartProd2 = np.array([posset for posset in CartProd if list(posset[0]) != list(posset[1])])

   #print(cartProd2)
   # Seperate into acetone atoms only array and MOF+Acetone array
   LAce = cartProd2[:, 1] # these are the Ace atoms
   LMOF = cartProd2[:, 0] # These are the MOF atoms

   #print(LAce)
   #print(LMOF[0])
   # distance between each atom and each acetone atom
   # Gets distances between each acetone and each MOF atom, as well as distances between each acetone atoms
   distCalc = lambda L1 : np.sqrt((L1[:, 0] - LMOF[:, 0]) ** 2 + (L1[:, 1] - LMOF[:, 1]) ** 2 + (L1[:, 2] - LMOF[:, 2]) ** 2)
   # return difference between calculated distance and desired distance to be minimized
   distc = distCalc(LAce)
   return np.array(distc) - dist

out = optimize.least_squares(fct, x0 = LocAce.flatten(), verbose = 1)
Ace_atoms = ['O', 'C', 'CH31', 'CH32']
print('Positions in the following order:',  *Ace_atoms)
print(*np.reshape(out['x'], (4,3)))