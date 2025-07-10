# 3D-AFM-codes
## The uploads of codes and raw data is still underway, if anything is wrong, please contact us.
This is a workflow to construct process the 3D-AFM raw data, do 3D-AFM averaging and plot the figures. The codes are developed in BIO-AFM-LAB at Weill Cornell Medicine under the supervision of Professor Simon Scheuring.

Developer: Runze Ma, Yining Jiang

Publication: Runze Ma, Cheng Giuseppe Chen, Yining Jiang, James Sturgis, Christophe Chipot and Simon Scheuring. Aquaporins destructure water above the pore. In submission.

User should email to the corresponding author of the paper for details about this work: Professor Simon Scheuring (sis2019@med.cornell.edu)

NOTE: Any usage of the codes should cite the publication mentioned above.

## System requirements:
1. Operating system for code development : Windows 10
2. Software for code development: MATLAB (MathWorks) 2024b, IgorPro 9.05 (optional, for 3D volume demonstration, e.g. supplementary video 1 and 2)
3. Additional add-ons: MIJI
4. Non-standard hardware: N/A

## Installation instructions: 
1. The codes require installation of MATLAB (MathWorks) 2024b. An installation guide can be found at: https://www.mathworks.com/help/install/.
2. IgorPro is recommended for 3D volume demonstration purpose. An installation guide can be found at: https://www.wavemetrics.com/.
3. MIJI is recommended (not required) for visualizing data. An installation guide can be found at: https://www.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab. If MIJI is not installed, user should comment out any code that uses MIJI for visualization (lines starting sith "MIJ.xxx").
4. The installation should take less than one hour on a "normal" desktop computer.

## General instructions:
These codes comprise of 3 parts, which are grouped into the subfolders of this repository. Instructions on using the codes are summarized here but can also be found in each subfolders.
### 1. raw_data_processing 
Note: These codes are developed to process the raw data saved from 3D-AFM setup.
#### main codes:
1. force_curve_3D_5_channels_stepZ_for_missing_lines.m
   Used for subtracting the trace scan data of 3D AFM, and calculated the corresponding force and energy curves.
2. force_curve_3D_5_channels_stepZ_for_missing_lines_retrace
   Used for subtracting the retrace scan data of 3D AFM, and calculated the corresponding force and energy curves.
3. display_3D_curve_all.mlapp
   Used for quick examine of the 3D matrices during and after the processing.
#### helper functions:
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
Still constructing...


### 3. plotting and demonstration
### figure plotting in main text
#### main codes:
1. figure_1_curve_calculation.m
	After loading the processed data, 1d spectra of amplitude, phase, energy, force will be given as plot at the given x-y pixel.
2. figure_2_image_rotate.m
	This is to rotate the 3D force data so that we can have a cross section that goes through the center of AqpZ in the 2d lattice.
3. figure_3_slices.m
	This involves rotating and cutting the 3D AFM average data to get slices 

### IgorPro 3D volume demonstration
#### main codes:
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