# Molecular Dynamics Trajectory Analysis

#### Scripts for analysis of LAMMPS trajectory outputs.

* Creating a heatmap of diffusing molecules in the simulation
	* COMDistribution.py - Creates a x,y,z coordinate file for the location of diffusing molecules, using the center of mass location for the molecule. 
	* diffusionHeatmap.py - Takes the output from COMDistribution.py and plots it as a heatmap using matplotlib. Plots the heatmap over a structure of the MOF (UiO-66)
	* MOFFrameworkTraj.py - Using a single snapshot of a molecule diffusing in a MOF, creates a seperate coordinate file of the MOF on its own. Loading this output in Ovito allows bonds to be drawn between molecules using pair-wise cutoffs. Rendering and saving an image of this serves as the background framework image for diffusionHeatmap.py.
* Calculating the fraction of Acetone molecules engaged in hydrogen bonding with mu3OHs in the MOF
	* hBondAnalysis.py - Calculates the fraction of diffusing molecules hydrogen bound to the MOF for any number of snapshots in a coordinate file
	* hydrogenBondAnalysis.py - Same as above, slightly different method.
