import pandas as pd
import numpy as np

# Starting locations on a set of mu3O-Zr-mu3OH
A = np.array([9.32453, 19.7188, 32.0558])
B = np.array([8.10115, 20.7263, 31.0265]) # Central atom
C = np.array([9.26509, 19.6265, 29.9703])







#A = np.array([21.7051, 21.6946, 21.7425]) #mu3O 2666

#B = np.array([22.6108, 20.4892, 20.5952]) # Zr, 2725

#C =	np.array([21.9727, 21.8991, 19.4837]) # mu3O, 2662

def calcAngle(A,B,C):
	ABVec = np.subtract(B,A)
	BCVec = np.subtract(C,B)

	dotProd = sum(ABVec * BCVec)

	MagAB = np.sqrt(sum(ABVec**2))
	MagBC = np.sqrt(sum(BCVec**2))

	angle = np.arccos([dotProd/(MagAB*MagBC)])*180/np.pi
	return 180 - angle

angleStart = calcAngle(A,B,C)
print(angleStart)

# Starting locations on a set of mu3O-Zr-mu3OH
A1 = np.array([29.8128, 9.12698, 21.9992])
B1 = np.array([31.0083, 10.4388, 22.6819]) # Central atom
C1 = np.array([30.8912, 10.5547, 20.739])

#4 (Omu3)	#
#A1 = np.array([29.8128, 9.12698, 21.9992])
#4 (Omu3)
#C1 = np.array([32.0773, 9.27783, 21.6769])
#5 (Zr)
#B1 = np.array([31.0083, 10.4388, 22.6819])

angleEnd = calcAngle(A1,B1,C1)
print(angleEnd)
