This folder contains all MATLAB code to do a drift correction on the 3D-AFM measurement based on lowest point mapping. Only run the main code and the main code must be run section by section.

## Main codes:
1. main_script_v1_AqpZ.m<br/>
	Based on the lowest point mapping(LPM) data, it calculates the drift constant and apply it to the energy matrix for later density construction.


## Helper functions:
2. drift_determination_v1.m<br/>
	Determine the drift constant based on the center of mass of the AqpZ.
3. particle_determination_v2b.m<br/>
	Use the drift information from code 1 to establish the x-y coordiantes relation before and after the drift correction.
4. energy_background_v3.m<br/>
	Use double exponential function to remove long-range energy background from the curves.


## Instructions:
1. Use the lowest point mapping(LPM) data to locate the center of mass of each AqpZ in the lattice. The data structure of the input is: lpmap_pp=[x1,y1;x2,y2;x3,y3;...], x-y coordiante of center of mass of each protein (in pixel), lpmap_np=[N1,N2;N2,N3;...],  neighboring situation of all the recorded center of mass in lpmap_pp, with N1,N2 being the neighboring centers, N starts from 0 in lpmap_pp. This step could be done by manually choosing the center or using cross-correlation method. An custom ImageJ plugin was used to locate the center of mass of AqpZ in the main text, which could be given upon reasonable request. An example coordinate was also given in the matlab codes to test the raw data.
2. For the section of "drift correction", initial guesses must be given for the calculation of drift constant.
3. For the section of "particle determination", default settings were given for the analysis of AqpZ. In these parameters, parsize is the radius of data with respect to the center of mass of each AqpZ, which would be considered into later calculation, res_tgt is the resolution of the sparse matrix for the data to be put into.
4. For the section of "energy background processing", zconst is provided as default for our raw data processing, which is the z piezo constant of the scanner.
5. Final output: a structure containing all the corrected data and correction information, prepared for the later processing.


## Raw data:
The original data can be found at: https://doi.org/10.5281/zenodo.15866205. For the main code here, the input data is "3.1 Input of 3D-AFM averaging.mat".
