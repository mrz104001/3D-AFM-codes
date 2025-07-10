This folder contains all the IgorPro scripts to display 3D matrix saved as tif/tiff file in 3D volumes. This helps plot the 3D volume demonstration in Figure 1 and Figure 2. You can also check the 3D data in other software if you don't have access to IgorPro software.

## Basic coding idea:
1. IgorPro provides a useful function called Gizmo to display 3D data as three orthogonal slices, easy to change color.
2. A simple panel is programmed in IgorPro to set the slice number in any of the three dimensions, change axis ratio, and generate a slice moving movie such as the supplementary movie 1 and 2.


## main codes:
1. GizmoSliceViewer.pxp
	Custom script base for the panels.
2. file_open.ipf
	Procedure file for the scripts.


## Instruction and reminder:
1. Open IgorPro, use File/Open Experiement to open GizmoSliceViewer.pxp file.
2. Use File/Open File/Procedure... to open file_open.ipf file, minimize or drag aside the file_open.ipf window.
3. Click 3D Plot/Gizmo for 3d afm 2, then click OK on the popped up Function Execution Error window.
4. Use the load button in the Gizmo 3D control panel to load a 3D data tif/tiff file and you can adjust the displayed slices.
5. A moving slice movie can also be generated via the movie button. More options are also given next to the button.
