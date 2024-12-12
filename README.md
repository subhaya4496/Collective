# Collective
File to execute-

ES2b_T.m -> Simulation with fixed/reflecting boundary conditions with and without wall interaction.

Parameter File -
Parameter_file.m -> Called in the executable files.

Boundary Conditions

Force Files -
1. ForceEB.m -> Called in the executable file. Calculates the total Elastic interactions between neighbouring cells. 
                It sums up all the neighbouring cell interactions.

 -  E_cell_updatedB.m -> Called inside ForceEB.m, where the cell interactions between 2 cells are calculated.

2. ForceLSB.m -> Called in the executable file. Calculates the linear spring/ steric interactions between overlapping cells.

3. wall_E.m -> Called in the executable file ES2b_T.m, calculates the boundary interactions (Steric interactions from the boundary is implemented here)

Verlet List-
LISTUB2.m -> Finds the neighbours within a fixed perimeter around the cell of interest. 
	        It is called inside ES2b_T.m, i.e. files with reflecting boundary conditions.
		(Requires some modification when it is reflecting parallel to one axis, but periodic parallel to the other.)

Defining Cell Position and Orientation-
diluteB.m -> Called in the executable files. Random or lattice arrangement of cells. Also controls the Box size or BoxL.

Saving Files for simulations-
savedataB.m -> Saves files in form *.mat and *.tif for data and images of simulations respectively.

Post-processing of images in ImageJ for movies.
