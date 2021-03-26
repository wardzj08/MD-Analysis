import matplotlib.pyplot as plt
import numpy as np

loadings = [0, 1, 2, 3, 5, 7]
hbond_frac = [0.534, 0.438, 0.48, 0.362, 0.256, 0.237]

plt.plot(loadings, hbond_frac)
plt.xlabel('Loading')
plt.ylabel('Average fraction of acetone hydrogen bound')
plt.savefig("hbondFracs.png", bbox_inches = 'tight')