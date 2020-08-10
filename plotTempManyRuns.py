import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


fig, ax = plt.subplots(1)
runs = [1, 2]
for runNum in runs:
	fRead = 'log.lammps'
	thermo = pd.read_csv(f'./{runNum}/{fRead}', skiprows = 1217, delim_whitespace=True,  nrows = 500001)
	print(thermo.head())
	print(thermo.tail())
	ax.plot(thermo['Step'], thermo['Temp'])
	ax.set_ylim(300, 350)


plt.savefig('thermoFig.png')

