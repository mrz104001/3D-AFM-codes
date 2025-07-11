This folder contains all MATLAB code for the plot of key panels in the main figure of the text. This is meant for the examination of all the functions we used for the data demonstration, so no detailed instruction will be provided here on how to use these codes. Only run the main code and the main code must be run section by section.

## Basic coding idea:
1. 3D matrix data are hard to be displayed in 2D figure, so MATLAB codes are used to plot 1d or 2d plots of the data.


## Main codes:
1. figure_1_curve_calculation.m
	After loading the processed data, 1d spectra of amplitude, phase, energy, force will be given as plot at the given x-y pixel.
2. figure_2_image_rotate.m
	This is to rotate the 3D force data so that we can have a cross section that goes through the center of AqpZ in the 2d lattice.
3. figure_3_slices.m
	This involves rotating and cutting the 3D AFM average data, to get vertical and horizontal slices.
4. Figure_4_MD_simulation.m
	This loads the MD simulation data, and cut the same cross section as in figure_3_slices.m
	
	
## Raw data:
The raw data of the 3D-AFM measurement on mica and AqpZ can be found at: https://doi.org/10.5281/zenodo.15856517. For main code 1, the input data set is 1.2. For main code 2, the input data is 2.2. For main code 3, the input data is 3.2. For main code 4, the input data could be 4.1-4.3.