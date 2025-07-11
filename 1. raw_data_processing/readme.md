This folder contains all the MATLAB code to process raw data from 3D-AFM acquired data via Zurich Intrument lock-in amplifier UHFLI.
The raw data is recorded as 4 data sets of amplitude(r), phase(phase), height(auxin1), X-scan signal(auxin2), a sample data on mica and on AqpZ could be found at the given URL: XXX.
Only run the main code and the main code must be run section by section!

## Basic coding idea:
1. Grouping the data into single approach and retract spectra on each x-y points of the 3D scan. The y cycles are determined via X-scan signal while the x cycles are determined by the height signal.
2. The grouped data is then saved into cell variables of amplitude, phase and height for later calculations, correction for the grouping and averaging for the cell variables are provided in the MATLAB code.
3. Then, the cell variable is aligned to the low height end of the spectra and transferred into 3D matrices of amplitude, phase and height.
4. Using Holsher's method, numerical integral and differential were deployed to get the force and energy matrices from amplitude, phase and height.


## Main codes:
1. force_curve_3D_5_channels_stepZ_for_missing_lines.m<br/>
   Used for subtracting the trace scan data of 3D AFM, and calculated the corresponding force and energy curves.
2. force_curve_3D_5_channels_stepZ_for_missing_lines_retrace.m<br/>
   Used for subtracting the retrace scan data of 3D AFM, and calculated the corresponding force and energy curves.
3. display_3D_curve_all.mlapp<br/>
   Used for quick examine of the 3D matrices during and after the processing.


## Helper functions:
4. double_expand.m<br/>
	Use double exponential to fit the amplitude and phase to larger height regions, used in force_cal_extended to calculate full range energy and force matrix.
5. find_period_max_fun_3D_cut.m<br/>
	Use sine wave fit to determine the x-scan period, used to help group the raw data into NxN cell variable.
6. force_cal_extended.m<br/>
	Used to calculate the full range energy and force matrix, since the amplitude and phase is expanded to longer range, data longer than the length of data obtained from force_cal.m is only considered for demonstration purpose.
7. force_cal.m<br/>
	Use Holsher method to calculate energy and force matrix.


## Instruction and reminder:
1. Grouping of the data might result in error of dividing between x-y point groups. For example, a supposed 50x50 x-y pixel 3D data, might give 49 or 51 data spectra in each x line while it is supposed to give 50. In this case, manual adjustment of the grouping is needed and the code is given in main code 1 (line 315-377) and 2 (line 313-372), which need to be used with caution.
2. Raw data from lock-in is given at URL: XXX, with the processed data file also available as reference and processing start for the next step.


## Raw data.
The raw data of the 3D-AFM measurement on mica and AqpZ can be found at: https://doi.org/10.5281/zenodo.15866205. 1.1 and 2.1 can be directly input into the main codes by running the loading section in the main code, 1.2 and 2.2 are the processed results.
