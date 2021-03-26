import matplotlib.pyplot as plt
import numpy as np

# Angle 9
K = 107.194822
C = [0.386699, -0.464206, 0.295566]

# Angle 14
K1 = 31.114298
K1 = 147.114298 # Modified
C1 = [0.286249, 0.225097, 0.262083]

theta = list(np.linspace(0,4*np.pi, 500))
doubletheta =  [2*x for x in theta]
E = K * (C[0] + C[1]*np.cos(theta) + C[2]*np.cos(doubletheta))
E1 = K1 * (C1[0] + C1[1]*np.cos(theta) + C1[2]*np.cos(doubletheta))

plt.plot(theta, E, label = 'Angle 9')
plt.plot(theta, E1, label = 'Angle 14')
plt.xlabel("Degrees (Radian)")
plt.ylabel("E (kcal/mol)")
plt.legend()
plt.savefig('tst2', bbox_inches = 'tight')
plt.show()
minE = np.min(E)
minE1 = np.min(E1)

locminE = [i for i, x in enumerate(E) if x <= minE+0.5]
locminE1 = [i for i, x in enumerate(E1) if x == minE1]
print(locminE)
print(locminE1)
minAnglesE = [theta[i] for i in locminE]
minAnglesE1 = [theta[i] for i in locminE1]

print(minAnglesE)
print(minAnglesE1)