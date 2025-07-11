#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later

// the procedure of loading 3d data, transfer it to the normal coordinate, do slices in the gizmo0 window.


// create the menu for the Gizmo 3d afm plot.
menu "3D Plot"
	"Gizmo for 3d afm 2",/Q, CreateGizmo3DPanel2()
end


Function CreateGizmo3DPanel2() : Panel
	killwindow/Z Gizmo3DPanel
	PauseUpdate; Silent 1		// building window...
	Dowindow/F Gizmo3DPanel
	variable w=500
	variable h=200
	variable i=1
	variable j=1
	variable base=0
	variable basew=0
	String/g Filename="Files to be loaded"
	String/g ColorName = "Turbo"
	variable/g XSliceN=0
	variable/g YSliceN=0
	variable/g ZSliceN=0
	variable/g Xmax=1
	variable/g Ymax=1
	variable/g Zmax=1
	variable/g Xsize
	variable/g Ysize
	variable/g Zsize
	variable/g method=1
	variable/g FPS=5, Direction=0, sleeptime=0.5
	variable/g FrameOnFlag = 1
	variable/g ColorInverseFlag = 1
	variable/g popNum = 1
	
	
	
	Display/K=1/W=(550,270,550+w,270+h)/N=Gizmo3DPanel as "Gizmo 3D control"
	
	// program the interface for the 3d plot.
	
	// make buttons.
	SetDrawEnv textyjust= 2,fname= "Segoe UI",fstyle= 1;DelayUpdate
	DrawText 10,10,"File being processed"
	
	Button ButtonLoadFile,pos={320,40},size={70,20},title="Load"
	Button ButtonLoadFile,fColor=(1,39321,39321) , proc=LoadFile2
	
	Button ButtonReLoad,pos={320,10},size={70,20},title="Reload"
	Button ButtonReLoad,fColor=(1,39321,39321), proc=ReloadButton
	
	
	TitleBox Filename, pos={10,35}, frame=5, fSize=12, fColor=(0,0,0), variable = $"Filename"
	TitleBox Filename fixedSize=1,size={300,25}
	
	//SetVariable SetAveraging, pos={20,70+i*30} ,size={140,20},title="Repeat Time" ,limits={1,1000,50},value=$"Averaging"
	
	
	//
	SetVariable ColorInverse, pos={10,10} ,size={110,20},title="Color Inverse?",limits={0,1,1},value=$"ColorInverseFlag"
	PopupMenu ColorName, pos={130,10} ,size={110,20},title="Color",value="Turbo;Grays;RedWhiteBlue",proc = popMenuColor
	
	// the axis tuning part.
	base=100
	i=0
	h=50
	GroupBox SliceControl title="Slice control",pos={10,base+i*h-25},size={380,175}
	
	Slider XSlice,pos={15,base+i*h},size={295,260},limits={0,9,1},vert=0,variable= $"XSliceN"
	SetVariable XSliceN, pos={315,base+i*h} ,size={70,20},title="X Slice",limits={1,10,0},value=$"XSliceN",proc=SliceNChange2
	i++
	Slider YSlice,pos={15,base+i*h},size={295,260},limits={0,9,1},vert=0,variable= $"YSliceN"
	SetVariable YSliceN, pos={315,base+i*h} ,size={70,20},title="Y Slice",limits={1,10,0},value=$"YSliceN",proc=SliceNChange2
	i++
	Slider ZSlice,pos={15,base+i*h},size={295,260},limits={0,9,1},vert=0,variable= $"ZSliceN"
	SetVariable ZSliceN, pos={315,base+i*h} ,size={70,20},title="Z Slice",limits={1,10,0},value=$"ZSliceN",proc=SliceNChange2
	i++
	
	
	// setting part, method of demonstrating the matrix, x y z ratio, etc.
	base=10
	i=0
	h=25
	GroupBox AspectRatio title="Aspect Ratio",pos={405,base+i*h},size={250,100}
	SetVariable Xsize, pos={415,base+i*h+25} ,size={100,20},title="X size(nm)",limits={0.1,200,0},value=$"Xsize"
	i++
	SetVariable Ysize, pos={415,base+i*h+25} ,size={100,20},title="Y size(nm)",limits={0.1,200,0},value=$"Ysize"
	i++
	SetVariable Zsize, pos={415,base+i*h+25} ,size={100,20},title="Z size(nm)",limits={0.1,10,0},value=$"Zsize"
	i++
	
	Button ButtonSetSize,pos={550,40},size={70,50},title="Set"
	Button ButtonSetSize,fColor=(1,39321,39321) , proc=SetSize2
	
	
	// create movie part.
	base=120
	i=0
	h=25
	GroupBox Movie title="Movie",pos={405,base+i*h},size={250,100}
	SetVariable FPS, pos={415,base+i*h+25} ,size={100,20},title="FPS(/s)",limits={1,20,0},value=$"FPS"
	i++
	SetVariable Direction, pos={415,base+i*h+25} ,size={100,20},title="Direction",help={"0 for y, 1 for x."},limits={0,1,1},value=$"Direction"
	i++
	SetVariable sleeptime, pos={415,base+i*h+25} ,size={100,20},title="sleep time(s)",limits={0.1,5,0},value=$"sleeptime"
	i++
	Button ButtonMovie,pos={550,150},size={70,50},title="Movie"
	Button ButtonMovie,fColor=(1,39321,39321) , proc=CreateMovie2
	
	
	// set the box frame part.
	Button ButtonFrame,pos={410,230},size={70,20},title="Frame"
	Button ButtonFrame,fColor=(1,39321,39321) , proc=FrameOn
	
	
	// set the paint.
	ModifyGizmo opName=background_black, operation=clearColor,data={0.0000,0.0000,0.0000,1.0000}
	ModifyGizmo opName=setsize, operation=scale,data={1,1,1}
	
	// reset the size of the matrix.
	variable xyz
	xyz=max(Xsize,Ysize,Zsize)
	ModifyGizmo opName=setsize, operation=scale,data={Xsize/xyz,Ysize/xyz,Zsize/xyz}
	
	// set the translation
	ModifyGizmo opName=translate, operation=translate,data={0,0,0}
	
	return 0
end


Function LoadFile2(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	switch( ba.eventCode )
		case 2: // mouse up
			SVAR Filename
			SVAR ColorName
			NVAR Xmax
			NVAR Ymax
			NVAR Zmax
			NVAR ZSliceN
			NVAR ColorInverseFlag
			wave waveAFM3d
		
			// wipe clean of the object list, the associated list would be removed as well.
			ModifyGizmo clearObjectList
			//ModifyGizmo clearDisplayList

			// loading 3d data. maybe add a button function later.
			ImageLoad/T=tiff/O/S=0/C=-1/Q/N='waveAFM3d'
			
			Xmax=DimSize('waveAFM3d',0);
			YMax=DimSize('waveAFM3d',1);
			ZMax=DimSize('waveAFM3d',2);
			ZSliceN=ZMax-1
			
			
			Slider XSlice,limits={0,Xmax-2,1}
			Slider YSlice,limits={0,Ymax-2,1}
			Slider ZSlice,limits={0,Zmax-1,1}
			
			Filename=S_fileName
			
			// plot it.
			ModifyGizmo stopUpdates
			AppendToGizmo surface=root:'waveAFM3d',name=ZPlane
			ModifyGizmo ModifyObject=ZPlane,objectType=surface,property={ srcMode,32}
			ModifyGizmo ModifyObject=ZPlane,objectType=surface,property={ plane,ZMax-1},property={ surfaceCTab,$ColorName}
			ModifyGizmo ModifyObject=ZPlane,objectType=surface,property={ inverseSurfaceCTAB,ColorInverseFlag},property={xStep, 2}, property={yStep, 2}
			AppendToGizmo surface=root:'waveAFM3d',name=YPlane
			ModifyGizmo ModifyObject=YPlane,objectType=surface,property={ srcMode,64},property={ surfaceCTab,$ColorName}
			ModifyGizmo ModifyObject=YPlane,objectType=surface,property={ inverseSurfaceCTAB,ColorInverseFlag},property={xStep, 2}, property={yStep, 2}
			AppendToGizmo surface=root:'waveAFM3d',name=XPlane
			ModifyGizmo ModifyObject=XPlane,objectType=surface,property={ srcMode,128},property={ surfaceCTab,$ColorName}
			ModifyGizmo ModifyObject=XPlane,objectType=surface,property={ inverseSurfaceCTAB,ColorInverseFlag},property={xStep, 2}, property={yStep, 2}
			AppendToGizmo Axes=BoxAxes,name=BoxFrame
			ModifyGizmo ModifyObject=BoxFrame,objectType=Axes,property={ -1,axisColor,1,1,1,1}
			ModifyGizmo ModifyObject=BoxFrame,objectType=Axes,property={ -1,lineWidth,1}
			
			// adding base.
			AppendToGizmo surface=root:'waveAFM3d',name=ZPlane0
			ModifyGizmo ModifyObject=ZPlane0,objectType=surface,property={ srcMode,32},property={ surfaceCTab,$ColorName}
			ModifyGizmo ModifyObject=ZPlane0,objectType=surface,property={ inverseSurfaceCTAB,ColorInverseFlag},property={xStep, 2}, property={yStep, 2}
			AppendToGizmo surface=root:'waveAFM3d',name=YPlane0
			ModifyGizmo ModifyObject=YPlane0,objectType=surface,property={ srcMode,64},property={ surfaceCTab,$ColorName}
			ModifyGizmo ModifyObject=YPlane0,objectType=surface,property={ inverseSurfaceCTAB,ColorInverseFlag},property={xStep, 2}, property={yStep, 2}
			ModifyGizmo ModifyObject=YPlane0,objectType=surface,property={ plane,YMax-1}
			AppendToGizmo surface=root:'waveAFM3d',name=XPlane0
			ModifyGizmo ModifyObject=XPlane0,objectType=surface,property={ srcMode,128},property={ surfaceCTab,$ColorName}
			ModifyGizmo ModifyObject=XPlane0,objectType=surface,property={ inverseSurfaceCTAB,ColorInverseFlag},property={xStep, 2}, property={yStep, 2}
			ModifyGizmo ModifyObject=XPlane0,objectType=surface,property={ plane,XMax-1}
			ModifyGizmo resumeUpdates

			// add display.
			ModifyGizmo stopUpdates
			ModifyGizmo setDisplayList=-1, object=ZPlane
			ModifyGizmo setDisplayList=-1, object=YPlane
			ModifyGizmo setDisplayList=-1, object=XPlane
			ModifyGizmo setDisplayList=-1, object=ZPlane0
			ModifyGizmo setDisplayList=-1, object=YPlane0
			ModifyGizmo setDisplayList=-1, object=XPlane0
			ModifyGizmo setDisplayList=-1, object=BoxFrame
			ModifyGizmo resumeUpdates
		
	endswitch

	return 0
end


Function SliceNChange2(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	switch( sva.eventCode )
		case 6: // value changed by dependency event.
			NVAR XSliceN
			NVAR YSliceN
			NVAR ZSliceN
			NVAR Xmax
			NVAR Ymax
			NVAR Zmax
			NVAR Xsize
			NVAR Ysize
			NVAR Zsize
			NVAR method
			wave waveAFM3d
			ModifyGizmo stopUpdates
			
			switch(method)
			case 0:	// default method, only plot one slice in each direction.
				ModifyGizmo ModifyObject=XPlane,objectType=surface,property={ plane,XSliceN}
				ModifyGizmo ModifyObject=YPlane,objectType=surface,property={ plane,YSliceN}
				ModifyGizmo ModifyObject=ZPlane,objectType=surface,property={ plane,ZSliceN}
				break
			case 1:	// maintaining a cube plot in the corner.
				duplicate/O/R = [XSliceN,Xmax-1][YSliceN,Ymax-1][0,ZSliceN] waveAFM3d waveAFM3d_cut
				ModifyGizmo stopUpdates
				ModifyGizmo ModifyObject=XPlane,objectType=surface,property={ srcWave,root:'waveAFM3d_cut'}
				ModifyGizmo ModifyObject=YPlane,objectType=surface,property={ srcWave,root:'waveAFM3d_cut'}
				ModifyGizmo ModifyObject=ZPlane,objectType=surface,property={ srcWave,root:'waveAFM3d_cut'}
				ModifyGizmo ModifyObject=ZPlane,objectType=surface,property={ plane,ZSliceN}
				ModifyGizmo resumeUpdates
				variable xyz
				xyz=max(Xsize,Ysize,Zsize)
				
				// keep the shape ratio the same.
				//ModifyGizmo opName=setsize, operation=scale,data={Xsize/xyz*(1-(XSliceN)/(Xmax-1)),Ysize/xyz*(1-YSliceN/Ymax),Zsize/xyz*(ZSliceN/Zmax)}
				
				
				// keep the original point the same.
				//ModifyGizmo opName=translate, operation=translate,data={-Xsize/xyz*XSliceN/(Xmax-1)/2,0,0}
				break
			endswitch
			
			ModifyGizmo resumeUpdates
		return 0
	endswitch
end

Function SetSize2(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	switch( ba.eventCode )
		case 2: // mouse up
			NVAR Xsize
			NVAR Ysize
			NVAR Zsize
			variable xyz
			xyz=max(Xsize,Ysize,Zsize)
			ModifyGizmo opName=setsize, operation=scale,data={Xsize/xyz,Ysize/xyz,Zsize/xyz}
		return 0
	endswitch
end


Function CreateMovie2(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	switch( ba.eventCode )
		case 2: // mouse up
			NVAR Ymax
			NVAR XSliceN
			NVAR YSliceN
			NVAR ZSliceN
			NVAR Xmax
			NVAR Ymax
			NVAR Zmax
			NVAR FPS
			NVAR Direction
			NVAR sleeptime
			wave waveAFM3d
			variable i
			variable Nmax
			dowindow/F Gizmo0	// need to bring the recorded window to the front for the newmovie recognization.
			NewMovie/F=(FPS);
			i=0
			do
				if (Direction==0)
					YSliceN=i
					Nmax=Ymax
				elseif (Direction == 1)
					XSliceN=i
					Nmax=Xmax
				endif
				duplicate/O/R = [XSliceN,Xmax-1][YSliceN,Ymax-1][0,ZSliceN] waveAFM3d waveAFM3d_cut
				ModifyGizmo stopUpdates
				ModifyGizmo ModifyObject=XPlane,objectType=surface,property={ srcWave,root:'waveAFM3d_cut'}
				ModifyGizmo ModifyObject=YPlane,objectType=surface,property={ srcWave,root:'waveAFM3d_cut'}
				ModifyGizmo ModifyObject=ZPlane,objectType=surface,property={ srcWave,root:'waveAFM3d_cut'}
				ModifyGizmo ModifyObject=ZPlane,objectType=surface,property={ plane,ZSliceN}
				ModifyGizmo resumeUpdates
				sleep/s sleeptime
				print i
				AddMovieFrame
				i++
				if (i>=Nmax)
					break
				endif
			while(i>-1)
			CloseMovie
		return 0
	endswitch
end

Function FrameOn(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	switch( ba.eventCode )
		case 2: // mouse up
		NVAR FrameOnFlag
		
		switch (FrameOnFlag)
			case 1:		// Frame box is already on.
				ModifyGizmo ModifyObject=BoxFrame,objectType=Axes,property={ -1,lineWidth,0}
				FrameOnFlag=0		// set it to be off
				return 0
			case 0:		// Frame box is already off.
				ModifyGizmo ModifyObject=BoxFrame,objectType=Axes,property={ -1,lineWidth,1}
				FrameOnFlag=1		// set it to be on.
				return 0
		endswitch
		
	endswitch
end


// popmenu action.
Function popMenuColor(s) : PopupMenuControl
    STRUCT WMPopupAction &s
   
    if (s.eventCode == 2)   // mouse up
    	SVAR ColorName
    	NVAR popNum
    	popNum = s.popNum
        ColorName = s.popStr
     endif
     return 0
end

Function ReloadButton(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	
	switch( ba.eventCode )
		case 2: // mouse up
			SVAR ColorName
			NVAR Xmax
			NVAR Ymax
			NVAR Zmax
			NVAR ZSliceN
			NVAR ColorInverseFlag
			wave waveAFM3d
			NVAR XSliceN
			NVAR YSliceN
			NVAR ZSliceN
			NVAR Xsize
			NVAR Ysize
			NVAR Zsize
		
			// plot it.
			ModifyGizmo stopUpdates
			ModifyGizmo ModifyObject=ZPlane,objectType=surface,property={ srcMode,32}
			ModifyGizmo ModifyObject=ZPlane,objectType=surface,property={ plane,ZMax-1},property={ surfaceCTab,$ColorName}
			ModifyGizmo ModifyObject=ZPlane,objectType=surface,property={ inverseSurfaceCTAB,ColorInverseFlag},property={xStep, 2}, property={yStep, 2}
			ModifyGizmo ModifyObject=YPlane,objectType=surface,property={ srcMode,64},property={ surfaceCTab,$ColorName}
			ModifyGizmo ModifyObject=YPlane,objectType=surface,property={ inverseSurfaceCTAB,ColorInverseFlag},property={xStep, 2}, property={yStep, 2}
			ModifyGizmo ModifyObject=XPlane,objectType=surface,property={ srcMode,128},property={ surfaceCTab,$ColorName}
			ModifyGizmo ModifyObject=XPlane,objectType=surface,property={ inverseSurfaceCTAB,ColorInverseFlag},property={xStep, 2}, property={yStep, 2}

			ModifyGizmo ModifyObject=BoxFrame,objectType=Axes,property={ -1,axisColor,1,1,1,1}
			ModifyGizmo ModifyObject=BoxFrame,objectType=Axes,property={ -1,lineWidth,1}
			
			// adding base.
			ModifyGizmo ModifyObject=ZPlane0,objectType=surface,property={ srcMode,32},property={ surfaceCTab,$ColorName}
			ModifyGizmo ModifyObject=ZPlane0,objectType=surface,property={ inverseSurfaceCTAB,ColorInverseFlag},property={xStep, 2}, property={yStep, 2}

			ModifyGizmo ModifyObject=YPlane0,objectType=surface,property={ srcMode,64},property={ surfaceCTab,$ColorName}
			ModifyGizmo ModifyObject=YPlane0,objectType=surface,property={ inverseSurfaceCTAB,ColorInverseFlag},property={xStep, 2}, property={yStep, 2}
			ModifyGizmo ModifyObject=YPlane0,objectType=surface,property={ plane,YMax-1}

			ModifyGizmo ModifyObject=XPlane0,objectType=surface,property={ srcMode,128},property={ surfaceCTab,$ColorName}
			ModifyGizmo ModifyObject=XPlane0,objectType=surface,property={ inverseSurfaceCTAB,ColorInverseFlag},property={xStep, 2}, property={yStep, 2}
			ModifyGizmo ModifyObject=XPlane0,objectType=surface,property={ plane,XMax-1}
			ModifyGizmo resumeUpdates

			// add display.
			ModifyGizmo stopUpdates
			ModifyGizmo setDisplayList=-1, object=ZPlane
			ModifyGizmo setDisplayList=-1, object=YPlane
			ModifyGizmo setDisplayList=-1, object=XPlane
			ModifyGizmo setDisplayList=-1, object=ZPlane0
			ModifyGizmo setDisplayList=-1, object=YPlane0
			ModifyGizmo setDisplayList=-1, object=XPlane0
			ModifyGizmo setDisplayList=-1, object=BoxFrame
			ModifyGizmo resumeUpdates
			
			
			// redraw the slices.
			duplicate/O/R = [XSliceN,Xmax-1][YSliceN,Ymax-1][0,ZSliceN] waveAFM3d waveAFM3d_cut
			ModifyGizmo stopUpdates
			ModifyGizmo ModifyObject=XPlane,objectType=surface,property={ srcWave,root:'waveAFM3d_cut'}
			ModifyGizmo ModifyObject=YPlane,objectType=surface,property={ srcWave,root:'waveAFM3d_cut'}
			ModifyGizmo ModifyObject=ZPlane,objectType=surface,property={ srcWave,root:'waveAFM3d_cut'}
			ModifyGizmo ModifyObject=ZPlane,objectType=surface,property={ plane,ZSliceN}
			ModifyGizmo resumeUpdates
			variable xyz
			xyz=max(Xsize,Ysize,Zsize)
		
	endswitch
end