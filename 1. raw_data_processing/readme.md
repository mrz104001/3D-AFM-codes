This folder contains all the MATLAB code to process raw data from 3D-AFM acquired data via Zurich Intrument lock-in amplifier UHFLI.
The raw data is recorded as 4 data sets of amplitude(r), phase(phase), height(auxin1), X-scan signal(auxin2), a sample data on mica and on AqpZ could be found at the given URL: XXX.

## Basic coding idea:
1. Grouping the data into single approach and retract spectra on each x-y points of the 3D scan. The y cycles are determined via X-scan signal while the x cycles are determined by the height signal.
2. The grouped data is then saved into cell variables of amplitude, phase and height for later calculations, correction for the grouping and averaging for the cell variables are provided in the MATLAB code.
3. Then, the cell variable is aligned to the low height end of the spectra and transferred into 3D matrices of amplitude, phase and height.
4. Using Holsher's method, numerical integral and differential were deployed to get the force and energy matrices from amplitude, phase and height.


## main codes:
1. force_curve_3D_5_channels_stepZ_for_missing_lines.m
   Used for subtracting the trace scan data of 3D AFM, and calculated the corresponding force and energy curves.
2. force_curve_3D_5_channels_stepZ_for_missing_lines_retrace
   Used for subtracting the retrace scan data of 3D AFM, and calculated the corresponding force and energy curves.
3. display_3D_curve_all.mlapp
   Used for quick examine of the 3D matrices during and after the processing.


## Instruction and reminder:
1. Grouping of the data might result in error of dividing between x-y point groups. For example, a supposed 50x50 x-y pixel 3D data, might give 49 or 51 data spectra in each x line while it is supposed to give 50. In this case, manual adjustment of the grouping is needed and the code is given in main code 1 and 2, which need to be used with caution.
2. Raw data from lock-in is given at URL: XXX, with the processed data file also available as reference and processing start for the next step.
