#pragma rtGlobals=1		// Use modern global access method.

Menu "UIUC"
	"Thermal Image Processing",ThermalipDriver()
End

Function ThermalipDriver()

	// If the panel is already created, just bring it to the front.
	DoWindow/F ThermalipPanel
	if (V_Flag != 0)
		return 0
	endif
	
	String dfSave = GetDataFolder(1)
	// Create a data folder in Packages to store globals.
	NewDataFolder/O/S root:packages:Thermalip
	
	//ReloadImageNames()
	Variable maxLimit = NumVarOrDefault(":gmaxLimit",NaN)
	Variable/G gmaxLimit= maxLimit
	Variable minLimit = NumVarOrDefault(":gminLimit",NaN)
	Variable/G gminLimit= minLimit
	Variable MaskDilations = NumVarOrDefault(":gMaskDilations",10)
	Variable/G gMaskDilations= MaskDilations
	Variable Threshold = NumVarOrDefault(":gThreshold",0.001)
	Variable/G gThreshold= Threshold	
	Variable LayerNum = NumVarOrDefault(":gLayerNum",0)
	Variable/G gLayerNum= LayerNum	
	Variable Transpose = NumVarOrDefault(":gTranspose",0)
	Variable/G gTranspose= Transpose	
	
	Variable Xstart = NumVarOrDefault(":gXstart",0)
	Variable/G gXstart= Xstart	
	Variable Xend = NumVarOrDefault(":gXend",1023)
	Variable/G gXend= Xend	
	Variable Ystart = NumVarOrDefault(":gYstart",0)
	Variable/G gYstart= Ystart	
	Variable Yend = NumVarOrDefault(":gYend",31)
	Variable/G gYend= Yend	
	
	Variable OffsetSignal = NumVarOrDefault(":gOffsetSignal",0)
	Variable/G gOffsetSignal = OffsetSignal	
	Variable scaleSignal = NumVarOrDefault(":gScaleSignal",1)
	Variable/G gScaleSignal = scaleSignal
	
	Variable layer2 = NumVarOrDefault(":gLayer2",0)
	Variable/G gLayer2 = layer2	
	Variable targLayer = NumVarOrDefault(":gTargLayer",1)
	Variable/G gTargLayer = targLayer	
	

	
	String /G gOperNames = "Arithmetic;Sensitivity;Replace data;Deconvolution"//;Beta"
	
	Make/O/N=(3) GRanges
	
	// Create the control panel.
	Execute "ThermalipPanel()"
	
	refreshElements(1)
	
	//Reset the datafolder to the root / previous folder
	SetDataFolder dfSave

End //SmartLithoDriver

Window ThermalipPanel(): Panel

	PauseUpdate; Silent 1		// building window...
	NewPanel /K=1 /W=(485,145, 780, 460) as "Thermal Image Processing"
	SetDrawLayer UserBack
	
	//Popupmenu imageselector,pos={17,72},size={135,18},title="Image", proc=ImageSelectorPM
	//Popupmenu imageselector,value= root:packages:Thermalip:Gimagenames,live= 1
	
	//Button button_RefreshNames,pos={196,110},size={150,20},title="Refresh Image List", proc=ReloadImageNames
	
	//Popupmenu layerselector,pos={17,12},size={135,18},title="Image", proc=LayerSelectorPM
	//Popupmenu layerselector,value= root:packages:Thermalip:Gimagenames,live= 1
		
	SetVariable SV_name,pos={20,15},size={255,18},title="Image Name", disable=2
	SetVariable SV_name,value= root:packages:MFP3D:Main:Display:LastTitle,live= 1
	
	Popupmenu pp_modType,pos={19,43},size={135,18},title="Operation",live= 1, proc=ProcOper	
	Popupmenu pp_modType,value= root:packages:Thermalip:gOperNames
	
	SetVariable SV_layernum,pos={178,48},size={98,18},title="Layer Index"
	SetVariable SV_layernum,value= root:packages:Thermalip:gLayerNum,live= 1
	
	Button button_undo,pos={20,73},size={256,23},title="Undo", proc=ThermalButtonFunc
	
	//////////////////////////////////////////////////////////////// ARITHMETIC ////////////////////////////////////////////////////////////////////////////
	
	//DrawText 18, 493, "Arithmetic Signal Modification:"
	
	SetVariable SV_AddVolts,pos={42,107},size={105,18},title="Add (mV)"
	SetVariable SV_AddVolts,value= root:packages:Thermalip:gOffsetSignal,live= 1
	
	Button button_Add,pos={173,104},size={75,20},title="Add", proc=addmV
	
	SetVariable SV_ScaleSignal,pos={42,137},size={105,18},title="Scale"
	SetVariable SV_ScaleSignal,value= root:packages:Thermalip:gScaleSignal,live= 1
	
	Button button_Scale,pos={173,137},size={75,20},title="Scale", proc=ScaleSignal
	
	Button button_invert,pos={45,174},size={85,30},title="Invert", proc=ThermalButtonFunc
	
	Button button_Flip,pos={160,174},size={85,30},title="Flip Vertical", proc=flipImage2
	
	SetVariable SV_maxlimit,pos={42,227},size={105,18},title="Max (mV)"
	SetVariable SV_maxlimit,value= root:packages:Thermalip:GMaxLimit,live= 1
	
	SetVariable SV_minlimit,pos={42,267},size={105,18},title="Min (mV)"
	SetVariable SV_minlimit,value= root:packages:Thermalip:GMinLimit,live= 1
	
	Button button_truncate,pos={174,224},size={75,57},title="Truncate", proc=ThermalButtonFunc
	
	//////////////////////////////////////////////////////////////// SENSITIVITY ////////////////////////////////////////////////////////////////////////////
	
	//DrawText 24, 249, "Manual Sensitivity Analysis:"
	
	ValDisplay VD_range1,pos={24,86},size={70,18}
	ValDisplay VD_range1,value= root:packages:Thermalip:Granges[0],live= 1
	
	ValDisplay VD_range2,pos={111,86},size={70,18}
	ValDisplay VD_range2,value= root:packages:Thermalip:Granges[1],live= 1
	
	ValDisplay VD_range3,pos={206,86},size={70,18}
	ValDisplay VD_range3,value= root:packages:Thermalip:Granges[2],live= 1
	
	Button Button_SetRange_0,pos={21,106},size={75,20},title="Set", proc=SetRangeFunc
	
	Button Button_SetRange_1,pos={109,106},size={75,20},title="Set", proc=SetRangeFunc
	
	Button Button_SetRange_2,pos={203,106},size={75,20},title="Set", proc=SetRangeFunc
		
	Button Button_GenerateString,pos={96,144},size={158,27},title="Generate!", proc=GetExcelString
	
	//SetVariable SV_maskdilations,pos={42,194},size={155,28},title="Mask Dilation Iterations"
	//SetVariable SV_maskdilations,value= root:packages:Thermalip:GMaskDilations,live= 1
	
	//Button button_stats,pos={212,191},size={91,21},title="Print Stats", proc=ThermalButtonFunc
	
	///////////////////////////////////////////////////////////// OVERWRITING /////////////////////////////////////////////////////////////////////////
	
	//DrawText 18, 335, "Single Layer Overwrite:"
	
	Checkbox CB_Transpose, pos = {42, 106}, size={10,10}, title="Transpose"// AFTER Truncating"
	Checkbox CB_Transpose, value=root:packages:Thermalip:gTranspose, live=1, proc=TransposeListener
	
	SetVariable SV_xstart,pos={45,136},size={98,18},title="X start"
	SetVariable SV_xstart,value= root:packages:Thermalip:gXstart,live= 1
	
	SetVariable SV_xend,pos={174,136},size={98,18},title="X end"
	SetVariable SV_xend,value= root:packages:Thermalip:gXend,live= 1
	
	SetVariable SV_Ystart,pos={45,165},size={98,18},title="Y start"
	SetVariable SV_Ystart,value= root:packages:Thermalip:gYstart,live= 1
	
	SetVariable SV_Yend,pos={174,165},size={98,18},title="Y end"
	SetVariable SV_Yend,value= root:packages:Thermalip:gYend,live= 1
	
	Button button_overWrite,pos={43,195},size={230,30},title="Load & Overwrite Layer", proc=OverwriteLayer
	
	
	
	////////////////////////////////////////////////////////// DECONVOLUTION /////////////////////////////////////////////////////////////////////////
	
	//DrawText 18, 335, "Single Layer Overwrite:"
		
	SetVariable SV_layer2,pos={45,113},size={98,18},title="Layer 2"
	SetVariable SV_layer2,value= root:packages:Thermalip:gLayer2,live= 1
	
	SetVariable SV_targetLayer,pos={174,113},size={98,18},title="Target Layer"
	SetVariable SV_targetLayer,value= root:packages:Thermalip:gTargLayer,live= 1
	
	Button button_Decon,pos={43,142},size={230,30},title="Layer 1 - Layer 2 = Target Layer", proc=DeconTopo
	
	
	//////////////////////////////////////////////////////////////// FILTERING ////////////////////////////////////////////////////////////////////////////
	
	//DrawText 16,117, "Edge Filter:"
	
	Button button_read,pos={40, 86},size={134,25},title="Read Filtered From Disk", proc=ReadImageFromDisk
	
	Button button_get1stDiff,pos={187, 86},size={117,25},title="Show 1st Differential", proc=ThermalButtonFunc
	
	SetVariable SV_threshold,pos={41,129},size={188,26},title="First Difference Threshold:"
	SetVariable SV_threshold,value= root:packages:Thermalip:GThreshold,live= 1
	
	Button button_filter, pos={242,125}, size={60,23}, title="Filter", proc=ThermalButtonFunc

		
	SetDrawEnv textrgb= (0,0,65280), fstyle= 1
	DrawText 122, 308, "Suhas Somnath, UIUC 2010"
	
EndMacro //ThermalipPanel

Function ProcOper(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	if( pa.eventCode == 2)
		 refreshElements(pa.popNum)	
	endif
	
End //ScanModeProc

Function refreshElements(mode)

	Variable mode;
		
		Variable isArith = mode==1;
		Variable isSens = mode==2;
		Variable isRepl = mode==3;
		Variable isFilt = mode==5;
		Variable isDecon = mode == 4;
		
		//Arithmetic;Sensitivity Analysis;Replace data;Beta"
		
		if(isSens)
			ModifyControl button_undo, disable= 1
		else
			ModifyControl button_undo, disable=0
		endif
		
		//ModifyControl SV_layernum, disable= !isRepl
		
		ModifyControl SV_AddVolts, disable= !isArith
		ModifyControl button_Add, disable= !isArith
		ModifyControl SV_ScaleSignal, disable= !isArith
		ModifyControl button_Scale, disable= !isArith		
		ModifyControl SV_maxlimit, disable= !isArith
		ModifyControl SV_minlimit, disable= !isArith
		ModifyControl button_truncate, disable= !isArith
		ModifyControl button_invert, disable= !isArith
		ModifyControl button_Flip, disable= !isArith
		
				
		ModifyControl SV_threshold, disable= !isFilt
		ModifyControl button_filter, disable= !isFilt
		ModifyControl button_read, disable= !isFilt
		ModifyControl button_get1stDiff, disable= !isFilt
						

		ModifyControl VD_range1, disable= !isSens
		ModifyControl VD_range2, disable= !isSens
		ModifyControl VD_range3, disable= !isSens
		ModifyControl Button_SetRange_0, disable= !isSens
		ModifyControl Button_SetRange_1, disable= !isSens
		ModifyControl Button_SetRange_2, disable= !isSens
		ModifyControl Button_GenerateString, disable= !isSens
		//ModifyControl SV_maskdilations, disable= !isSens
		//ModifyControl button_stats, disable= !isSens
		

		ModifyControl CB_Transpose, disable= !isRepl
		ModifyControl SV_xstart, disable= !isRepl
		ModifyControl SV_xend, disable= !isRepl
		ModifyControl SV_Ystart, disable= !isRepl
		ModifyControl SV_Yend, disable= !isRepl
		ModifyControl button_overWrite, disable= !isRepl
		
		ModifyControl SV_layer2, disable = !isDecon;
		ModifyControl SV_targetLayer, disable = !isDecon;
		ModifyControl button_Decon, disable = !isDecon;
		
End

Function SetRangeFunc(ctrlname) : ButtonControl
	String ctrlname
	
	Variable index = -1
	
	Variable RemInd = FindLast(CtrlName,"_")
	if (RemInd > -1)
		CtrlName = CtrlName[RemInd+1,Strlen(CtrlName)-1]
		index = str2num(CtrlName)
	else
		print "Error in Button function"
		return -1
	endif
	
	String dfSave = GetDataFolder(1)
	
	SetDataFolder root:packages:Thermalip
	
	Wave GRanges
	
	SetDataFolder root:packages:MFP3D:Main:Analyze:Section
	
	Wave SectionWaveI
	
	GRanges[index] = 1000*abs(SectionWaveI[0] - SectionWaveI[1]) // getting them in mV not V
	
	SetDataFolder dfSave
End

Function GetExcelString(ctrlname) : ButtonControl
	String ctrlname
	
	String dfSave = GetDataFolder(1)
	
	SetDataFolder root:packages:Thermalip
	
	Wave GRanges
	
	String output = "=Average("+num2str(GRanges[0])+","+num2str(GRanges[1])+","+num2str(GRanges[2])+")"
	
	Granges[0] = 0
	Granges[1] = 0
	Granges[2] = 0
		
	print output
	
	SetDataFolder dfSave
	
End


Function ThermalButtonFunc(ctrlname) : ButtonControl
	String ctrlname
	
	Variable RemInd = FindLast(CtrlName,"_")
	if (RemInd > -1)
		CtrlName = CtrlName[RemInd+1,Strlen(CtrlName)-1]
	else
		print "Error in Button function"
		return -1
	endif
	
	String dfSave = GetDataFolder(1)
	
	SetDataFolder root:packages:MFP3D:Main:Display
	
	SVAR LastTitle
	Variable index = strsearch(LastTitle, " ", 0)
	if(index < 0)
		DoAlert 0, "No such Image!"
		return 0;
	endif
	
	String imgname = LastTitle[0,index-1]
	Variable layernum = -1
	
	String GraphName = StringFromList(0,WinList(cOfflineBaseName+"*",";","WIN:1"))	//get the name of the top graph
	if (strlen(GraphName) == 0)		//anything there?
		DoAlert 0, "No such Image!"
		return 0							//nope
	endif
	
	String DataFolder, ImageName
	GetGraphData(GraphName,DataFolder,ImageName,LayerNum)
	
	
	SetDataFolder root:packages:Thermalip
	
	NVAR GMaxLimit, GMinLimit, GThreshold, GMaskDilations
			
	if(validateImage(imgname,layernum))
		return -1
	endif
	
	strswitch (ctrlName)

		case "invert":
			copyFromImage(imgname, layernum, "ImgBackup")
			invertgraph(imgname,layernum)
		break
		
		case "truncate":
			copyFromImage(imgname, layernum, "ImgBackup")
			truncategraph(imgname,layernum, GMaxLimit*1e-3, GMinLimit*1e-3)
		break
		
		case "stats":
			AutoImageStats(GMaskDilations)
		break
		
		case "filter":
			copyFromImage(imgname, layernum, "ImgBackup")
			filterImage(GThreshold)
			writeToImage(imgname, layernum,"Filtered")
		break
		
		case "undo":
			writeToImage(imgname, layernum,"ImgBackup")
		break
		
		case "get1stDiff":
			DoFlattenFunc("DoFlatten_0")
			copyFromImage(imgname, layernum, "ImgBackup")
			makeFirstDiff()
			showDiff(0)
		break
		
	endswitch
	
	LocalDisplayAutoRangeFunc(GraphName)
	
	SetDataFolder dfSave
	
End

function validateImage(wavenam,layer)
	String wavenam
	Variable layer	
	
	SetDataFolder root:Images
	
	Wave image = $wavenam
	
	if(!WaveExists(image))
		DoAlert 0, "Wave not found\nCheck Name of Wave please!"
		SetDataFolder dfSave
		return -1
	endif
	
	if(layer > DimSize(image3,2) || layer < 0)
		DoAlert 0, "Layer not found\nTip: Layer indices start from 0 not 1"
		SetDataFolder dfSave
		return -1
	endif
	
	return 0	
end

function writeToImage(imgname, layer,sourcewave)

	String imgname
	Variable layer
	String sourcewave
	
	String dfSave = GetDataFolder(1)
	SetDataFolder root:Images
	
	Wave image = $imgname
	
	SetDataFolder root:packages:Thermalip
	
	Wave inwave = $sourcewave
	if(!exists(sourcewave))
		print "Destination wave not recognized"
		return -1;
	endif
	
	Variable i , j
	for (i=0; i < DimSize(image,0); i += 1)
		for (j=0; j< DimSize(image,1); j += 1)
			image[i][j][layer] = inwave[i][j]
		endfor
	endfor
	
	SetDataFolder dfSave

End

function copyFromImage(imgname, layer, destwave)

	String imgname
	Variable layer
	String destwave
	
	String dfSave = GetDataFolder(1)
	SetDataFolder root:Images
	
	Wave image = $imgname
	
	SetDataFolder root:packages:Thermalip
	
	if(!exists(destwave))
		print "Error: no such wave already existing. Making a new one"
		Make/O/N=(DimSize(image,0),DimSize(image,1)) $destwave
	else
		Redimension/N=(0,0) $destwave
		Make/O/N=(DimSize(image,0),DimSize(image,1)) $destwave
	endif
	
	Wave outwave = $destwave
	
	Variable i , j
	for (i=0; i < DimSize(image,0); i += 1)
		for (j=0; j< DimSize(image,1); j += 1)
			outwave[i][j] = image[i][j][layer]
		endfor
	endfor
	
	SetDataFolder dfSave

End

function makeFirstDiff()

	String dfSave = GetDataFolder(1)
	SetDataFolder root:packages:Thermalip
	
	Wave ImgBackup
	
	Make/O/N=(DimSize(ImgBackup,0),DimSize(ImgBackup,1)) Diff1
    
       Variable i=0
       Variable j=0
       
   	for(j=0; j<DimSize(ImgBackup,1); j+=1)
      		Diff1[0][j] = 0;
      		for(i=1; i<DimSize(ImgBackup,0);i+=1)
            		Diff1[i][j] = ImgBackup[i][j]-ImgBackup[i-1][j];
        	endfor
   	endfor
   	
   	SetDataFolder dfSave

end

function showDiff(colnum)
	Variable colnum
	
	String dfSave = GetDataFolder(1)
	SetDataFolder root:packages:Thermalip
	
	If(exists("diff1"))
		if(DimSize(ImgBackup,1) < colnum || colnum < 0)
			print "Error: Invalid column number"
			return -1
		endif
		Wave Diff1
		Display/K=1 Diff1[][colnum]
	endif
	
	SetDataFolder dfSave
	
End

function filterImage(threshold)
	Variable threshold;
	
	String dfSave = GetDataFolder(1)
	SetDataFolder root:packages:Thermalip
	
	makeFirstDiff()
	
	Wave Diff1
	
	Variable i=0
       Variable j=0
	
	for(j=0; j<DimSize(Diff1,1); j+=1)
	
		Variable startindex = -1;
		Variable total = 0;
		
		do
			//print "Looking at point i = " + num2str(i) + ", j = " + num2str(j) +" value = " + num2str(diff1[i][j])
			// Case 1 - start of peak:
            		if(abs(diff1[i][j]) > abs(threshold))
                		if(startindex < 0)
                			//print "Found start index at " + num2str(i) + ", " + num2str(j);
                    		startindex = i;
                    		total = diff1[i][j]
                		else
                    		total = total + diff1[i][j]
                		endif
            		endif
            		
            		// Case 2 - End of peak:
           		if(abs(diff1[i][j]) < abs(threshold) && startindex > 0)
                
               		Variable endindex = -1;
                		// Case a - single point - need to find missing ones:
				
				if(startindex == i-1)
                    		// Should start looking at the next 2-3 points
                    		// an opposite value should be waiting
                    		Variable last = min(i+3,DimSize(Diff1,0));
                    		Variable k = i
                    		for(k=i;k<last; k+=1)
                            		total = total + diff1[k][j]
                        			if(abs( diff1[k][j]) > abs(threshold))
                            			endindex = k;
                            			i = k;
                            			//print "Found new end index for single point at " + num2str(i) + ", " + num2str(j);
                            			break; 
                        			endif
                    		endfor
                    		//print "endindex = " + num2str(endindex)
                    		//return -1
                		else
                		
                			// Case b - multiple points found:
                    		// Try looking ahead two points?
                    		
                    		//print "Case 2b end point at "  + num2str(i-1) + ", " + num2str(j);
                    		endindex = i-1;
                    		
                		endif
                		
                		// Actual Filtering:
                		if(endindex > 0)
                		
                			Variable pos
                		
	                    	if(total > 0)
      		                  		// Rising
             	           		pos = endindex;
                   			else
                        			pos = startindex;
                   			endif
                      
                    		if(pos > startindex)
                    		
                    			//print "Setting points i = "  + num2str(startindex) + " to  " + num2str(pos-1) + ", j = "+ num2str(j) + ") = 0";                    		
                        			
                        			for (k=startindex; k< pos; k+=1)
                           			diff1[k][j] = 0;
                        			endfor
                        			     
                   			endif

					//print "Setting peak (i = "  + num2str(pos) + ", j = "+ num2str(j) + ") = " + num2str(total) ;
                   			diff1[pos][j] = total;
                   			
                   			                   	
                   			if(pos+1 < endindex)     
                   				//print "Setting points i = "  + num2str(pos+1) + " to  " + num2str(endindex) + ", j = "+ num2str(j) + ") = 0";         			
                    			for(k=pos+1; k<=endindex; k+=1)
                       				diff1[k][j] = 0; 
                    			endfor
                    		endif
                    		//reset start index
                			startindex = -1;
                			total = 0;
                		
                		endif

            		endif
            		
            		i+= 1; 

		while( i < DimSize(Diff1,0))
		
		i = 0;
		
		//print "Next j----------------------------------------------------------------------------------"
			
	endfor
	
	// Calculating the actual data:
    	Make/O/N=(DimSize(Diff1,0),DimSize(Diff1,1)) Filtered
    	
    	Wave ImgBackup
    	
   	for(j=0; j< DimSize(Diff1,1); j+=1)
      		Filtered[0][j] = ImgBackup[0][j]
      		for(i=1; i< DimSize(Diff1,0); i+=1)
          		Filtered[i][j] =  Filtered[i-1][j] + diff1[i][j];
       	endfor
    	endfor

	SetDataFolder dfSave

end

//The objective of this function is 
function innerFilter(threshold, tolmult)
	Variable threshold, tolmult
	
	String dfSave = GetDataFolder(1)
	SetDataFolder root:packages:Thermalip
	
	Wave Diff1
	
	Variable i=0
       Variable j=0
	
	for(j=0; j<DimSize(Diff1,1); j+=1)
		Variable startindex = -1;
		Variable total = 0;
		
		do
			// Case 1 - peak:
            		if(abs(diff1[i][j]) > abs(threshold)*tolmult)
                		if(startindex == -1)
                			startindex = i
                		endif
                	else
                	// Case 2 - not peak
                		// Case 2a - End of peak
                		if(startindex != -1)
						
                		endif
            		endif
            		
      
			
		
		
		while( i < DimSize(Diff1,0))
	endfor
	SetDataFolder dfSave
	
End

function truncategraph(wavenam,layer, maxval, minval)

	String wavenam
	Variable layer, maxval,minval	

	SetDataFolder root:Images
	
	Wave image = $wavenam
	WaveStats/Q image
	
	if(maxval == NaN)
		maxval = V_max+1;
	endif
	if(minval == NaN)
		minval = V_min-1
	endif
	
	Variable i , j
	for (i=0; i < DimSize(image,0); i += 1)
		for (j=0; j< DimSize(image,1); j += 1)
			if(image[i][j][layer] > maxval)
				image[i][j][layer] = maxval
			endif
			if(image[i][j][layer] < minval)
				image[i][j][layer] = minval
			endif
		endfor
	endfor
	
end

function MakeVcant(wavenam,Vtotlayer, Vsenselayer)

	String wavenam
	Variable Vtotlayer, Vsenselayer 	

	SetDataFolder root:Images
	
	Wave image = $wavenam
	WaveStats/Q image
		
	Variable i , j
	for (i=0; i < DimSize(image,0); i += 1)
		for (j=0; j< DimSize(image,1); j += 1)
			image[i][j][Vtotlayer] -= image[i][j][Vsenselayer]
		endfor
	endfor
	
end

function MakeVcantVcont(wavenam, Vtotalvalue, Vsenselayer)

	String wavenam
	Variable Vtotalvalue, Vsenselayer 	

	SetDataFolder root:Images
	
	Wave image = $wavenam
	WaveStats/Q image
		
	Variable i , j
	for (i=0; i < DimSize(image,0); i += 1)
		for (j=0; j< DimSize(image,1); j += 1)
			image[i][j][Vsenselayer] = Vtotalvalue - image[i][j][Vsenselayer]
		endfor
	endfor
	
end

function invertgraph(wavenam,layer)

	String wavenam
	Variable layer	

	SetDataFolder root:Images
	
	Wave image = $wavenam
	WaveStats/Q image
	
	Variable i , j
	for (i=0; i < DimSize(image,0); i += 1)
		for (j=0; j< DimSize(image,1); j += 1)
			image[i][j][layer] = V_max - image[i][j][layer]
			
		endfor
	endfor
	
end

Function ReloadImageNames()

	String dfSave = GetDataFolder(1)
	SetDataFolder root:packages:Thermalip
	SVAR GImageNames
	Variable i
	
	Wave namewave = root:Images:MeMListWave
	
	if(numpnts(namewave) > 0)
		//Execute("print namewave[0]")
		GImageNames = num2str(namewave[0])
			
		for(i=1;i<numpnts(namewave);i+=1)
			if(namewave[i] != 0)
				String mytemp=""
				Execute("print namewave[i]")
				print mytemp
				GImageNames += ";" + num2str(namewave[i])
			endif 							
		endfor	
	else
		GImageNames = ""
	endif
	
	SetDataFolder dfSave

End


Function ReadImageFromDisk(ctrlname): ButtonControl
	String ctrlname
	String oldSaveFolder = GetDataFolder(1)
	setdatafolder root:packages:Thermalip
	Variable refNum
	String outputPath
	Open /R /Z=2 /M="Select the text file containing the Filtered Image" refNum as ""
	if(refNum == 0)
		print "No file was open!"
		//return -1
	endif
	if (V_flag == -1)
		Print "Open cancelled by user."
		return -1
	endif
	if (V_flag != 0)
		DoAlert 0, "Error Opening file"
		return V_flag
	endif
	outputPath = S_fileName
	
	
	//Lets extract the pattern's name from its filename
	String filename = outputpath
	Variable enditer = FindLast(filename,":")
	filename = filename[enditer+1,strlen(filename)-1]
	//removing the .txt
	filename = RemoveEnding(filename , ".TXT")	
	//We cant allow any spaces within the wave's name:
	filename = ReplaceString(" ",filename,"_",1)
	
	if(cmpstr(outputPath,"")==0)
		DoAlert 0, "\t\tError!!\n\n\tYou did not choose any file!"
		return -1
	else
		//print "Wave name = " + filename
		readWaves(refNum,filename)
	endif
	setdatafolder oldSaveFolder
End //LoadWavesFromDisk

Function readWaves(filePointer,name)
	Variable filePointer
	String name
	
	// storing the waves directly into the LithoWaves
	String oldSaveFolder = GetDataFolder(1)
	SetDataFolder root:packages:ThermalIP
	
	String str
	//First line has the size of the waves:
	FReadLine filePointer, str
	Variable numrows, numcols, i, j
	
	sscanf str, "%d\t%d", numrows, numcols
		
	//print "# of lines" + num2str(wavesize) + ">> X = " + num2str(firstx) + ", Y = " + num2str(firsty)
	
	Make/O /N=(numrows,numcols) Filtered
		
	// Providing a handle now
	//Wave Xwave = $("X"+name)
		
	//Now looping over the remaining size to start filling in the wave:
	for(j=0; j<numcols; j= j+1)
		for (i=0; i<numrows; i= i+1)
			FReadLine filePointer, str
			Variable temp
			sscanf str, "%f", temp
			Filtered[i][j] = temp
		endfor
	endfor
	
	SetDataFolder oldSaveFolder
	
	Close filePointer
	return 0
End //End of ReadWaves

Function AutoImageStats(maskDilations)
	Variable maskDilations
	DoRoughnessFunc("DoRoughness_0")
	Wave imgdetails = root:Packages:MFP3D:Main:Analyze:Roughness:FullStats
	Variable range = imgdetails[%V_Max] -  imgdetails[%V_Min]
	range = range * 1000
	OfflineCheckFunc("FillMaskBox_3",1)
	AROfflineButtonFunc("ResetMask_3")
	OfflineCheckFunc("InverseMaskBox_3",0)
	Variable uppernoise = 1000*getRoughness(maskDilations)
	OfflineCheckFunc("InverseMaskBox_3",1)
	AROfflineButtonFunc("ResetMask_3")
	Variable lowernoise = 1000*getRoughness(maskDilations)
	OfflineCheckFunc("InverseMaskBox_3",0)
	AROfflineButtonFunc("ResetMask_3")
	//Variable avgnoise = (uppernoise + lowernoise)/2
	String outpt = num2str(uppernoise) + "\t" + num2str(lowernoise) + "\t" + num2str(range)
	//DoAlert 0, "Image stats in mV:\n" + outpt
	print "Upper noise,       Lower Noise,       Range"
	print outpt
End

Function getRoughness(maskDilations)
	Variable maskDilations
	OfflineImagePopFunc("MaskCalcTypePopup_3",2, "iterative")
	AROfflineButtonFunc("MakeMask_3")
	Variable i=0
	for(i=0; i<maskDilations; i+=1)	
		AROfflineButtonFunc("DilateMask_3")
	endfor
	DoRoughnessFunc("DoRoughness_0")
	Wave maskedimg = root:Packages:MFP3D:Main:Analyze:Roughness:MaskStats
	return maskedimg[%V_Sdev]
End

Function TransposeListener(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba
	
	String dfSave = GetDataFolder(1)
	SetDataFolder  root:packages:ThermalIP
	NVAR gTranspose

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			gTranspose = checked
			break
	endswitch
	
	SetDataFolder dfSave

	return 0
End

Function OverwriteLayer(ctrlname) : ButtonControl
	String ctrlname
	
	killWaves/Z FiltData0, TruncData
	
	String oldSaveFolder = GetDataFolder(1)
	SetDataFolder root:packages:ThermalIP
	NVAR gLayerNum, gTranspose, gXstart, gXend, gYstart, gYend
	
	// 1. Load data from disk
	LoadWave/Q/J/M/D/A=FiltData/K=0 //"C:Users:somnath2:Desktop:test.txt"
	if(exists("FiltData0")==0)
		print "Error: Did not load data. Aborting"
		return -1;
	endif
	Wave FiltData0
	
	// 2. Truncate data
	//(9,1088-1)
	//Duplicate/O/R=(gXstart,gXend)(gYstart,gYend) FiltData0 TruncData
	Duplicate/O/R=[gXstart,gXend] FiltData0 TruncData	
	// 3. Acquire target wave:
	SVAR imgraw = root:packages:MFP3D:Main:Display:LastTitle;
	Variable endindx =  strsearch(imgraw,"#0",0)
	if(endindx < 0)
		print "Error:No image name found. Try selecting the image. "
		return -1
	endif
	
	imgraw = imgraw[0,endindx-2];
	SetDataFolder root:images
	Wave targetWave = $(imgraw);
	
	// 4. Overwrite data!
	
	
	Variable i,j;
	
	for(i=0;i< DimSize(targetWave, 0);i=i+1)
		for(j=0; j< DimSize(targetWave, 1);j=j+1)
			if(gTranspose==1)
			 	targetWave[i][j][gLayerNum] = TruncData[j][i];
			 else
			 	 targetWave[i][j][gLayerNum] = TruncData[i][j];
			 endif
		endfor
	endfor
	
	// If you don't kill the waves here, it retains the old values for some reason.
	KillWaves TruncData, FiltData0;
	
	SetDataFolder oldSaveFolder
End

Function ScaleSignal(ctrlname) : ButtonControl
	String ctrlname
	
	String oldFolder = GetDataFolder(1)
	
	SetDataFolder root:packages:ThermalIP
	NVAR gLayerNum, gScaleSignal;
	
	if(gScaleSignal == 1)
		print "Nothing to scale"
		return 0;
	endif
	
	SVAR imgraw = root:packages:MFP3D:Main:Display:LastTitle;
	Variable endindx =  strsearch(imgraw,"#0",0)
	if(endindx < 0)
		print "Error:No image name found. Try selecting the image. "
		return -1
	endif
	
	imgraw = imgraw[0,endindx-2];
	SetDataFolder root:images
	Wave targetWave = $(imgraw);
	
	// 2. Scale data!
	
	Variable i,j;
	
	for(i=0;i< DimSize(targetWave, 0);i=i+1)
		for(j=0; j< DimSize(targetWave, 1);j=j+1)
			 targetWave[i][j][gLayerNum] = targetWave[i][j][gLayerNum]*gScaleSignal;
		endfor
	endfor
	
	SetDataFolder oldFolder
	
End

Function addmV(ctrlname) : ButtonControl
	String ctrlname
	
	String oldFolder = GetDataFolder(1)
	
	SetDataFolder root:packages:ThermalIP
	NVAR gLayerNum, gOffsetSignal;
	
	if(gOffsetSignal == 0)
		print "Nothing to Offset"
		return 0;
	endif
	
	SVAR imgraw = root:packages:MFP3D:Main:Display:LastTitle;
	Variable endindx =  strsearch(imgraw,"#0",0)
	if(endindx < 0)
		DoAlert 0,"Error:No image name found. Try selecting the image. "
		return -1
	endif
	
	imgraw = imgraw[0,endindx-2];
	SetDataFolder root:images
	Wave targetWave = $(imgraw);
	
	// 2. Offset data!
	
	Variable i,j;
	
	for(i=0;i< DimSize(targetWave, 0);i=i+1)
		for(j=0; j< DimSize(targetWave, 1);j=j+1)
			 targetWave[i][j][gLayerNum] = targetWave[i][j][gLayerNum] + gOffsetSignal*0.001;
		endfor
	endfor
	
	SetDataFolder oldFolder
	
End

Function flipImage2(ctrlname) : ButtonControl
	String ctrlname
	
	String oldFolder = GetDataFolder(1)
	
	SetDataFolder root:packages:ThermalIP
	NVAR gLayerNum;
	
	SVAR imgraw = root:packages:MFP3D:Main:Display:LastTitle;
	Variable endindx =  strsearch(imgraw,"#0",0)
	if(endindx < 0)
		DoAlert 0,"Error:No image name found. Try selecting the image. "
		return -1
	endif
	
	imgraw = imgraw[0,endindx-2];
	SetDataFolder root:images
	Wave targetWave = $(imgraw);
	Duplicate/O targetWave tempWave	
	// 2. Flip Image
	// [pixel][line][layer]
	// /R=[gXstart,gXend]
	
	Variable pix,lin;
	
	for(lin=0;lin< DimSize(targetWave, 1);lin=lin+1)
		for(pix=0;pix< DimSize(targetWave, 0);pix=pix+1)
			 targetWave[pix][DimSize(targetWave, 1)-1-lin][gLayerNum] = tempWave[pix][lin][gLayerNum];
		endfor
	endfor
	
	killWaves/Z tempWave
	
	SetDataFolder oldFolder
	
End

Function DeconTopo(ctrlname) : ButtonControl
	String ctrlname
	
	String oldFolder = GetDataFolder(1)
	
	SetDataFolder root:packages:ThermalIP
	NVAR gLayerNum,gLayer2,gTargLayer;	
	
	if(gLayerNum == gLayer2 || gLayer2 == gTargLayer)
		DoAlert 0, "Pick different valid source and/or target layer numbers"
		SetDataFolder oldFolder
		return -1;
	endif
	
	SVAR imgraw = root:packages:MFP3D:Main:Display:LastTitle;
	Variable endindx =  strsearch(imgraw,"#0",0)
	if(endindx < 0)
		DoAlert 0,"Error:No image name found. Try selecting the image. "
		SetDataFolder oldFolder
		return -1
	endif
	
	imgraw = imgraw[0,endindx-2];
	SetDataFolder root:images
	Wave targetWave = $(imgraw);
	
	if(gLayerNum >= DimSize(targetWave,2) || gLayerNum >= DimSize(targetWave,2) || gLayerNum >= DimSize(targetWave,2))
		SetDataFolder oldFolder
		return -1;
	endif
	
	// 2. Offset data!
	
	Variable i,j;
	
	for(i=0;i< DimSize(targetWave, 0);i=i+1)
		for(j=0; j< DimSize(targetWave, 1);j=j+1)
			 targetWave[i][j][gTargLayer] = targetWave[i][j][gLayerNum] - targetWave[i][j][gLayer2];
		endfor
	endfor
	
	SetDataFolder oldFolder
	
End