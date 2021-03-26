# "Heatmap" plotter for guest molecule position data. Makes use of a kernal density estimator to plot a continuous function
# of the distribution of adsorbate molecules based on their center of masses. Displays the distribution over a still image of the MOF.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.image as mpimg

df = pd.read_csv('dump.COMNH3_1216.lammpstrj', skiprows = 9, names=['id', 'mol', 'type', 'x', 'y', 'z'], delim_whitespace=True)
sns.set(rc={'figure.figsize':(20,20)})
framework_img = mpimg.imread('./data/1900Frame1.png') # read in still framework file
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
# Control heaatmap: vary n_levels, gridsize and bw to best display data. BW is the smoothing factor, increasing it makes the density plot more homogenous
# as it approaches 0, individual molecules in the system stand out from one another
htmp = sns.kdeplot(df.y, df.z,
                 cmap="Reds", shade=True, shade_lowest= False, n_levels = 25, gridsize = 100, bw = 4, ax = ax)#, cut=0) # bw=.15,
#im.collections[0].set_alpha(0)
print(max(df.z))
print(max(df.y))
plt.imshow(framework_img, aspect = htmp.get_aspect(), zorder =1,
          extent = [0, 64, 0, 47])#im.get_xlim() + im.get_ylim(), zorder = 1)
# For no labels keep following two lines in
plt.xlabel("")
plt.ylabel("")
# for no tick marks, keep following two lines in
plt.xticks([])
plt.yticks([])
plt.show()
plt.savefig('UiO67_NH3Ads2_P1216_Smooth4.png', bbox_inches = 'tight')
#plt.pcolormesh(X,Y,Z)




