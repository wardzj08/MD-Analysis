import pandas as pd
import numpy as np

# 4 atoms in the dihedral
P = [1, 0, 1] # A1
Q = [0, 0, 0] # A2
R = [1, 0, 0] # A3
S = [0, 1, 1] # A4

P = [19.5391, 11.5115, 29.8893] #2666 (mu3-O) (collapsed)
Px = [21.8617, 11.5115, 32.2119] #2666 (mu3-O) (uncollapsed)
Q = [20.7004, 12.832, 31.0506] #2725 #zr
R = [19.5391, 11.5115, 32.2119] #2658 #mu3oh
S = [18.2186, 10.3502, 31.0506] #2716 #zr (trans)
Sx = [20.7004, 10.3502, 33.5324] #2716 #zr (cis)

# 2 planes are A1-A2-A3 and A4-A3-A2

# Good
def calc_Plane(P, Q, R):
	PQVec = np.subtract(P,Q)
	PRVec = np.subtract(P,R)

	cp = np.cross(PQVec, PRVec)
	d = -(cp[0] * P[0] + cp[1] * P[1] + cp[2] * P[2])
	#print(d)
	#print(cp)
	plane_coords = list(cp) + [d]
	return plane_coords

#coords1 = calc_Plane(P, Q, R)
#coords2 = calc_Plane(S, R, Q)
#print(coords1)
#print(coords2)

#numerator = coords1[0]*coords2[0] + coords1[1]*coords2[1] + coords1[2]*coords2[2]
#denominator = np.sqrt(coords1[0]**2 + coords1[1]**2 + coords1[2]**2) * np.sqrt(coords2[0]**2 + coords1[1]**2 + coords2[2]**2)

#dihedral = np.arccos(numerator/denominator) * 180/np.pi
#print(dihedral)



#b1 = np.subtract(P, Q)
#b2 = np.subtract(Q, R)
#b3 = np.subtract(R,S)

#numer = np.dot(np.cross(b1, b2), np.cross(b3, b2))
#denom = np.abs(np.cross(b1, b2)) * np.abs(np.cross(b3, b2))

#angle1 = np.arccos(numer/denom)
#print(angle1)



def wiki_dihedral(p):
    """formula from Wikipedia article on "Dihedral angle"; formula was removed
    from the most recent version of article (no idea why, the article is a
    mess at the moment) but the formula can be found in at this permalink to
    an old version of the article:
    https://en.wikipedia.org/w/index.php?title=Dihedral_angle&oldid=689165217#Angle_between_three_vectors
    uses 1 sqrt, 3 cross products"""
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    b0xb1 = np.cross(b0, b1)
    b1xb2 = np.cross(b2, b1)

    b0xb1_x_b1xb2 = np.cross(b0xb1, b1xb2)

    y = np.dot(b0xb1_x_b1xb2, b1)*(1.0/np.linalg.norm(b1))
    x = np.dot(b0xb1, b1xb2)

    return np.degrees(np.arctan2(y, x))

angle1 = wiki_dihedral(np.array([P, Q, R, S]))
print(angle1)

P1 = [20.7328, 10.0595, 30.9964] #2666 (mu3-O) (collapse)
P1x = [21.8039, 11.2134, 32.1712] #2666 (mu3-O) (uncollapse)
Q1 = [21.0698, 9.98895, 29.0553] #2725 #zr
R1 = [20.5899, 10.0505, 31.1803] #2658 #mu3oh
S1 = [18.6244, 10.272, 31.1899] #2716 #zr
S1 = [20.5188, 10.2213, 33.1421] #2716 #zr
angle2 = wiki_dihedral(np.array([P1, Q1, R1, S1]))
print(angle2)