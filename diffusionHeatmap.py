# "Heatmap" plotter from diffusing molecule position data. Makes use of a kernal density estimator to plot a continuous function
# of the distribution of adsorbate molecules based on their center of masses. Displays the distribution over a still image of the MOF.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.image as mpimg

df = pd.read_csv('dump.COMAce1Load.lammpstrj', skiprows = 9, names=['id', 'mol', 'type', 'x', 'y', 'z'], delim_whitespace=True)
#print(df.iloc[1])

sns.set(rc={'figure.figsize':(10,10)})
framework_img = mpimg.imread('UiO66frame2Nosim.png')
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
htmp = sns.kdeplot(df.x, df.y,
                 cmap="Reds", shade=True, shade_lowest= False, n_levels = 100, gridsize = 100, ax = ax) # bw=.15,
#im.collections[0].set_alpha(0)
plt.imshow(framework_img, aspect = htmp.get_aspect(),
          extent = [0, 41.4008, 0, 41.4008], zorder = 1)#im.get_xlim() + im.get_ylim(), zorder = 1)
# For no labels keep following two lines in
plt.xlabel("")
plt.ylabel("")
# for no tick marks, keep following two lines in
plt.xticks([])
plt.yticks([])
plt.show()
plt.savefig('MoleculeDistributionAce1LoadinUiO66TESTnosimnoTicks.png', bbox_inches = 'tight')
