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
#dist = [1.8, 2.94, 4.27, 4.25, 4.37]
AceInterDist = [1.229, 2.40154392, 2.40154392, 1.229, 1.51999949, 1.51999949, 2.40154392, 1.51999949, 2.594794, 2.40154392, 1.51999949, 2.594794]
dist = dist + AceInterDist
#dist = np.zeros(3)
# O C CH31 CH32
# Initial xyz values for acetone atoms (will help keep in line distance between these atoms)
xBead = [0, 0] #[0, 0, 0, 0]
yBead = [1.297397, -1.297397]#[0, 0, 1.297397, -1.297397]
zBead = [2.020934, 2.020934]#]][0, 1.229, 2.020934, 2.020934]
#xBead = [23.30983974, 24.31376747]#, 24.24619224, 25.69395881]
#yBead = [17.63103319, 17.18174926]#, 16.28809239, 17.49837908]
#zBead = [23.80153849, 24.3431314]#, 25.56672571, 23.81048819]

# 4 cluster atoms from the center Zr cluster in the super cell
# H(atom # 2581) Zr1 (atom # 2717) Zr2 (atom # 3633) Zr3 (atom # 3625)
xMOF = [23.00983974, 24.03583974]#[21.96126, 20.70040, 20.70040, 23.18217]
yMOF = [18.25103319, 17.65103319]#[18.75353, 18.21863, 20.70040, 20.70040]
zMOF = [24.00153849, 24.3103849]#[22.63941, 20.70040, 23.18217, 20.70040]
dist = [2.40154392, 2.40154392, 1.51999949, 1.51999949, 2.594794, 2.594794]

# Put atoms into matrix
LocAce = np.transpose(np.array((xBead, yBead, zBead))) # Ace
LocMOF = np.transpose(np.array((xMOF, yMOF, zMOF))) # MOF

# Function to minimize (difference of two distances)
def fct(LOCAce):
	# reshape to xyz coordinates for each atom instead of a list of all coordinates
   LOCAce = np.reshape(LOCAce, (2,3)) # Ace

   # array of all acetone and MOF reference atoms
   LOCS = np.concatenate((LocMOF, LOCAce))
   # calculate cartesian product of Acetone atoms and MOF atoms + Acetone atoms
   CartProd = np.array(list(itertools.product(*[LOCS, LOCAce])))
   # Remove coordinate sets that are two positions for the same atom (ex acetone oxygen paired with acetone oxygen)
   cartProd2 = np.array([posset for posset in CartProd if list(posset[0]) != list(posset[1])])
  # cartProd2 = CartProd
   print("hello")
   print(cartProd2)
   # Seperate into acetone atoms only array and MOF+Acetone array
   LAce = cartProd2[:, 1] # these are the Ace atoms
   LMOF = cartProd2[:, 0] # These are the MOF atoms + Ace atoms
   #print(LAce[1], LMOF[1])
   # distance between each atom and each acetone atom
   # find distances between each acetone and each MOF atom, as well as distances between each acetone atoms
   distCalc = lambda L1 : np.sqrt((L1[:, 0] - LMOF[:, 0]) ** 2 + (L1[:, 1] - LMOF[:, 1]) ** 2 + (L1[:, 2] - LMOF[:, 2]) ** 2)
   # return difference between calculated distance and desired distance to be minimized
   #print(distCalc(LAce) - dist)
   #print(dist[1])
   x =  distCalc(LAce) - dist
   objective = (x**2).sum()
   return objective

# least squares optimization of the distance function. Minimize the difference between the calculated distances and the desired distances
#out = optimize.least_squares(fct, x0 = LocAce.flatten(), verbose = 1)#, #xtol = 1e-25, ftol = 1e-25, gtol = 1e-25)#, #method = 'lm')
out = optimize.minimize(fct, x0 = LocAce.flatten(), method = "SLSQP")
Ace_atoms = ['O', 'C', 'CH31', 'CH32']
print('Positions in the following order:',  *Ace_atoms)
print(*np.reshape(out['x'], (2,3)))

