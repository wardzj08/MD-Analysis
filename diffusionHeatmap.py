# "Heatmap" plotter from diffusing molecule position data. Makes use of a kernal density estimator to plot a continuous function
# of the distribution of adsorbate molecules. Displays the distribution over an image of the MOF/framework. 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.image as mpimg 

df = pd.read_csv('dump.COM.lammpstrj', skiprows = 9, names=['id', 'mol', 'type', 'x', 'y', 'z'], delim_whitespace=True)
print(df.iloc[1])
# Compress into 2 dimensions by adding x,y; y,z; x,z
#totaled = pd.DataFrame({'x': df.x, 'y': df.y})
#totaled = totaled.append(pd.DataFrame({'x': df.y, 'y': df.z}), ignore_index=True)
#totaled = totaled.append(pd.DataFrame({'x': df.x, 'y': df.z}), ignore_index=True)
#print(len(totaled))
sns.set(rc={'figure.figsize':(10,10)})
framework_img = mpimg.imread('UiO66frame1.png')
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
htmp = sns.kdeplot(df.x, df.y,
                 cmap="Reds", shade=True, shade_lowest= False, n_levels = 100, gridsize = 100, ax = ax) # bw=.15,
#im.collections[0].set_alpha(0)
plt.imshow(framework_img, aspect = htmp.get_aspect(),
          extent = [0, 41.4008, 0, 41.4008], zorder = 1)#im.get_xlim() + im.get_ylim(), zorder = 1)
plt.show()
plt.savefig('MoleculeDistributioninUiO66.png', bbox_inches = 'tight')
