# Plotting functions for WindowMeasure.py
# Given the atom pair names, plots the change is distance with respect to time for each atom pair on an individual subplot
# The colors list allows the colors on each plot to be customized to match colors of the atoms in a vizualization of the transition state.
# Thus a side by side vizual of these plots with the transition state can show how each atom pair set moves as the adsorbate moves through the window.
# While the trend is fairly visible for each atom set when using all the data, a smoothing filter has been applied to reduce the noise on this data in an attempt
# to make the trends painfully obvious.

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

# read in the distance file
distances12 = np.loadtxt('BoydNoAceDistancesAveragedInOut12.csv', delimiter = ',')
distances13 = np.loadtxt('BoydNoAceStartDistancesAveragedInOut13.csv', delimiter = ',')
distances23 = np.loadtxt('BoydNoAceStartDistancesAveragedInOut23.csv', delimiter = ',')
#print(distances)

# Plot titles and the corresponding colors

#plts = ['1C', '2C', '3C', '4C', '5C', '6C']
#titles = ['1C', '2C', '3C', '4C', '5C', '6C']
#colors = ['#ffaaff', '#00aaff', '#ffaa7f', '#55aa7f', '#aaaaff', '#fceea7']

plts = ['1-2 In', '1-2 Out', '1-3 In', '1-3 Out', '2-3 In', '2-3 Out']

# Plotting function.
# For each atom pair, plots the distance v. time on an individual subplot
def plotPlts(plts):
    fig = plt.figure()
    fig.set_size_inches(18,14)
    for i,df in enumerate(plts):
    	# Select subplot
        ax = fig.add_subplot(2,3, i+1)

        # Plot distances with a smoothing filter applied
        ax.plot(savgol_filter(distances[:, i], 5, 2), color = colors[i])
        ax.set_title(titles[i])

    return fig

fig = plt.figure()
fig.set_size_inches(18,14)
ax = fig.add_subplot(1,3, 1)
ax.plot(distances12[:, 0], label = '1-2 In')
ax.plot(distances12[:, 1], label = '1-2 Out')
ax.legend(loc = "upper right", fontsize = 20)
ax.set_title("1-2 Linkers", fontsize = 24)
ax = fig.add_subplot(1, 3, 2)
ax.plot(distances13[:, 0], label = '1-3 In')
ax.plot(distances13[:, 1], label = '1-3 Out')
ax.legend(loc = "upper right", fontsize = 20)
ax.set_title("1-3 Linkers", fontsize = 24)
ax = fig.add_subplot(1, 3, 3)
ax.plot(distances23[:, 0], label = '2-3 In')
ax.plot(distances23[:, 1], label = '2-3 Out')
ax.legend(loc = "upper right", fontsize = 20)
ax.set_title("2-3 Linkers", fontsize = 24)
#final = plotPlts(plts)
fig.savefig('BoydNoAceWindowDistAveragedPairs2.png')