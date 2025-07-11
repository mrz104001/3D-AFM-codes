# 3D-AFM-codes
## Primitive version of the code and raw data has been uploaded. Please contact us if anything is unclear.
This is a workflow to process the 3D-AFM raw data, do 3D-AFM averaging and plot the figures. The codes are developed in BIO-AFM-LAB at Weill Cornell Medicine under the supervision of Professor Simon Scheuring.

Developer: Runze Ma, Yining Jiang

Publication: Runze Ma, Cheng Giuseppe Chen, Yining Jiang, James Sturgis, Christophe Chipot and Simon Scheuring. Aquaporins destructure water above the pore. In submission.

User should email to the corresponding author of the paper for details about this work: Professor Simon Scheuring (sis2019@med.cornell.edu)

NOTE: Any usage of the codes should cite the publication mentioned above.

## System requirements:
1. Operating system for code development : Windows 10
2. Software for code development: MATLAB (MathWorks) 2024b, ImageJ, IgorPro 9.05 (optional, for 3D volume demonstration, e.g. supplementary video 1 and 2)
3. Additional add-ons: MIJI
4. Non-standard hardware: N/A

## Installation instructions: 
1. The codes require installation of MATLAB (MathWorks) 2024b. An installation guide can be found at: https://www.mathworks.com/help/install/.
2. IgorPro is recommended for 3D volume demonstration purpose. An installation guide can be found at: https://www.wavemetrics.com/.
3. MIJI is recommended (not required) for visualizing data. An installation guide can be found at: https://www.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab. If MIJI is not installed, user should comment out any code that uses MIJI for visualization (lines starting sith "MIJ.xxx").
4. The installation should take less than one hour on a "normal" desktop computer.

## General instructions:
These codes comprise of 3 parts, which are grouped into the subfolders of this repository.
Instructions on using the codes are summarized here but can also be found in each subfolders.
### 1. Raw_data_processing 
Note: These codes are developed to process the raw data saved from 3D-AFM setup.
#### Main codes:
1. force_curve_3D_5_channels_stepZ_for_missing_lines.m
   Used for subtracting the trace scan data of 3D AFM, and calculated the corresponding force and energy curves.
2. force_curve_3D_5_channels_stepZ_for_missing_lines_retrace
   Used for subtracting the retrace scan data of 3D AFM, and calculated the corresponding force and energy curves.
3. display_3D_curve_all.mlapp
   Used for quick examine of the 3D matrices during and after the processing.
#### Helper functions:
4. double_expand.m
	Use double exponential to fit the amplitude and phase to larger height regions, used in force_cal_extended to calculate full range energy and force matrix.
5. find_period_max_fun_3D_cut.m
	Use sine wave fit to determine the x-scan period, used to help group the raw data into NxN cell variable.
6. force_cal_extended.m
	Used to calculate the full range energy and force matrix, since the amplitude and phase is expanded to longer range, data longer than the length of data obtained from force_cal.m is only considered for demonstration purpose.
7. force_cal.m
	Use Holsher method to calculate energy and force matrix.
#### Instruction and reminder:
1. Grouping of the data might result in error of dividing between x-y point groups. For example, a supposed 50x50 x-y pixel 3D data, might give 49 or 51 data spectra in each x line while it is supposed to give 50. In this case, manual adjustment of the grouping is needed and the code is given in main code 1 (line 315-377) and 2 (line 313-372), which need to be used with caution.
2. Raw data from lock-in is given at URL: XXX, with the processed data file also available as reference and processing start for the next step.


### 2. 3D AFM averaging
### 2.1 pre-processing
#### Main codes:
1. main_script_v1_AqpZ.m
	Based on the lowest point mapping(LPM) data, it calculates the drift constant and apply it to the energy matrix for later density construction.
#### Helper functions:
1. drift_determination_v1.m
	Determine the drift constant based on the center of mass of the AqpZ.
2. particle_determination_v2b.m
	Use the drift information from code 1 to establish the x-y coordiantes relation before and after the drift correction.
3. energy_background_v3.m
	Use double exponential function to remove long-range energy background from the curves.
## Instructions:
1. Use the lowest point mapping(LPM) data to locate the center of mass of each AqpZ in the lattice. The data structure of the input is: lpmap_pp=[x1,y1;x2,y2;x3,y3;...], x-y coordiante of center of mass of each protein (in pixel), lpmap_np=[N1,N2;N2,N3;...],  neighboring situation of all the recorded center of mass in lpmap_pp, with N1,N2 being the neighboring centers, N starts from 0 in lpmap_pp. This step could be done by manually choosing the center or using cross-correlation method. An custom ImageJ plugin was used to locate the center of mass of AqpZ in the main text, which could be given upon reasonable request. An example coordinate was also given in the matlab codes to test the raw data.
2. For the section of "drift correction", initial guesses must be given for the calculation of drift constant.
3. For the section of "particle determination", default settings were given for the analysis of AqpZ. In these parameters, parsize is the radius of data with respect to the center of mass of each AqpZ, which would be considered into later calculation, res_tgt is the resolution of the sparse matrix for the data to be put into.
4. For the section of "energy background processing", zconst is provided as default for our raw data processing, which is the z piezo constant of the scanner.
5. Final output: a structure containing all the corrected data and correction information, prepared for the later processing.

### 2.2 Density construction
#### Main codes:
1. volume_construction_essemble_v3_AqpZ.m
	Merge the drift-corrected energy data above multiple AqpZ into one single dataset, and apply a kernel generated by the helper function to connect the isolated points into a surface.
2. force_calculation_v1.m
	Use numerical differential to calculate the force from energy matrix.
3. force_adjustment_v1.m
	Restore the values of the matrix to real value since the kernel application step actually decrease the value by spreading the value of the point to the volume of the kernel.
## Helper function.
1. make_3D_LAFM_kernel1a.m
	Generate a kernel to be used in the main code 1.
## Instructions.
1. Run "volume_construction_essemble_v3_AqpZ.m" first, it needs the input of the structure data set from pre-processing. The parameters of the first section is the default used in the main text.
2. The result of "volume_construction_essemble_v3_AqpZ.m" will be saved at the end of the script, the filename of which needs to be input into "force_calculation_v1.m" for the continue processing.
3. Lastly, run "force_adjustment_v1.m", no extra input is needed. "ppvolumeF3" will be the final result of 3D-AFM averaged data.


### 3. Plotting and demonstration
### 3.1 Figure plotting in main text
#### Main codes:
1. figure_1_curve_calculation.m
	After loading the processed data, 1d spectra of amplitude, phase, energy, force will be given as plot at the given x-y pixel.
2. figure_2_image_rotate.m
	This is to rotate the 3D force data so that we can have a cross section that goes through the center of AqpZ in the 2d lattice.
3. figure_3_slices.m
	This involves rotating and cutting the 3D AFM average data to get slices 

### 3.2 IgorPro 3D volume demonstration
#### Main codes:
1. GizmoSliceViewer.pxp
	Custom script base for the panels.
2. file_open.ipf
	Procedure file for the scripts.
#### Instruction and reminder:
1. Open IgorPro, use File/Open Experiement to open GizmoSliceViewer.pxp file.
2. Use File/Open File/Procedure... to open file_open.ipf file, minimize or drag aside the file_open.ipf window.
3. Click 3D Plot/Gizmo for 3d afm 2, then click OK on the popped up Function Execution Error window.
4. Use the load button in the Gizmo 3D control panel to load a 3D data tif/tiff file and you can adjust the displayed slices.
5. A moving slice movie can also be generated via the movie button. More options are also given next to the button.


## Raw data
The raw data for all the codes here are deposited here: https://doi.org/10.5281/zenodo.15856517. Please refer to the readme file in each subfolder for detailed description.