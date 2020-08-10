import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
img = plt.imread("rigidTST2.png")
fig, ax = plt.subplots()
df = pd.read_csv('FinalTetraEnergyProfile.csv')
ax.imshow(img)
#ax.spines['top'].set_visible(False)
#ax.spines['left'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
#ax.spines['right'].set_visible(False)
#ax.set_xticks([])
#ax.set_yticks([])
#xy2imgxy = lambda x,y: (img.size[0] * x / np.max(ticklx), img.size[1] * (np.max(tickly) - y) / np.max(tickly))
#ticklx = np.linspace(0,10)
#tickly = np.linspace(0,17)
#tickpx,tickpy = xy2imgxy(ticklx,tickly)
ax.imshow(img, aspect='auto')

plt.savefig('FINAL.png')
import seaborn as sns
import matplotlib.image as mpimg


sns.set(rc={'figure.figsize':(10,10)})
framework_img = mpimg.imread('rigidTST2.png')
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
pt = ax.plot(df['#Coor'], df['Free'], zorder = 1)
ax.set_xlim(-1, 11)
#im.collections[0].set_alpha(0)

ax.imshow(img, zorder=0, extent=[0, 10.0, 0.0, 17.0])
# For no labels
plt.xlabel("")
#plt.ylabel("")
# for no tick marks
plt.xticks([])
#plt.yticks([])
plt.show()
plt.savefig('FINAL1.png', bbox_inches = 'tight')
