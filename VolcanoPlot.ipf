#pragma TextEncoding = "MacRoman"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <WaveSelectorWidget>
#include <PopupWaveSelector>

// Figure in the TPD54 paper was generated with
// https://github.com/quantixed/TPD54/commit/accb6a86619aa4668e5dadac655e55630e9a55f2
// Load data into Igor and then run the command from the Macro menu.
// Note that naming of the columns is important

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////
Menu "Macros"
	SubMenu "Proteomics"
		"Load MaxQuant Data...", /Q, LoadMaxQuantData()
		"Volcano Plot...", /Q, VolcanoIO_Panel()
		"PCA Only...", /Q, MakePCAWaveSelectorPanel()
		"Label Top 10", /Q, LabelTopTenWorkflow()
	end
End

////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////

Function LoadMaxQuantData()
	if(LoadMaxQuantFile() == 0)
		VolcanoIO_Panel()
	endif
End

Function VolcanoWorkflowWrapper(STRING prefix1,STRING prefix2,VARIABLE baseVal,VARIABLE pairOpt)
	MakeVolcano(prefix1,prefix2,baseVal,pairOpt)
	MakeColorTableWave()
	MakeVPlot()
	TableInterestingValues()
	AddSignificantHitsToVolcano()
	MakeMeanComparison()
	FromVolcanoToPCA()
	MakeTheLayout()
End

Function LabelTopTenWorkflow()
	LabelTopXProts(10)
End

////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////

Function LoadMaxQuantFile()
	NewDataFolder/O/S root:data 
	LoadWave/A/D/J/K=0/L={0,0,0,0,0}/W/O/Q ""
	if (strlen(S_Path) == 0) // user pressed cancel
		return -1
	endif
	WAVE/Z/T Gene_names, Protein_names
	Duplicate/O Gene_names, root:SHORTNAME
	Duplicate/O Protein_names, root:NAME
	String wList = WaveList("LFQ_Intensity*",";","")
	Variable nWaves = ItemsInList(wList)
	String wName
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		Duplicate/O $wName, $("root:" + ReplaceString("LFQ_Intensity_",wName,""))
	endfor
	
	SetDataFolder root:
	return 0
End

// This function drives the whole program
/// @param	prefix1	string prefix that will select all waves of condition1
/// @param	prefix2	string prefix that will select all waves of condition1
/// @param	baseVal	proteins that are absent have this value
/// @param  pairOpt	1 is paired, 0 is not
Function MakeVolcano(prefix1,prefix2,baseVal,pairOpt)
	String prefix1,prefix2
	Variable baseVal,pairOpt
	
	String wList1 = WaveList(prefix1,";","")
	String wList2 = WaveList(prefix2,";","")
	if(ItemsInList(wList1) == 0 || ItemsInList(wList2) == 0)
		DoAlert 0, "Missing data"
		return -1
	endif
	// make sure the lists are in order
	wList1 = SortList(wList1)
	wList2 = SortList(wList2)
	// Possibly we should use a more sophisticated way to check the groups match?
	Concatenate/O wList1, allCond1
	Concatenate/O wList2, allCond2
	
	// deal with baseVal
	TransformImputeBaseVal(allCond1,baseVal)
	TransformImputeBaseVal(allCond2,baseVal)
	
	// now do T-tests
	Variable nProt = dimsize(allCond1,0)
	Make/O/N=(nProt) allTWave,colorWave=0
	Variable pVar
	
	Variable i
	
	for(i = 0; i < nProt; i += 1)
		MatrixOp/O/FREE w0 = row(allCond1,i) ^ t
		MatrixOp/O/FREE w1 = row(allCond2,i) ^ t
		if(pairOpt == 0)
			StatsTTest/Q/Z w0,w1
		else
			StatsTTest/Q/Z/PAIR w0,w1
		endif
		WAVE/Z W_StatsTTest
		if(V_flag == 0)
			if(pairOpt == 0)
				pVar = W_StatsTTest[9] // p-value
			else
				pVar = W_StatsTTest[6] // p-value for Paired
			endif
		else
			pVar = 1
		endif
		allTWave[i] = pVar
	endfor
	
	// make mean waves - these need transformation back
	allCond1[][] = 10^(allCond1[p][q])
	allCond2[][] = 10^(allCond2[p][q])
	MatrixOp/O meanCond1 = sumrows(allCond1)
	meanCond1 /=ItemsInList(wList1)
	MatrixOp/O meanCond2 = sumrows(allCond2)
	meanCond2 /=ItemsInList(wList2)
	
	// ratio wave
	MatrixOp/O ratioWave = meanCond1 / meanCond2
	
	Duplicate/O ratioWave,ratioWave_log2
	ratioWave_log2 = log(abs(ratioWave[p])) / log(2)
	// assign colors
	colorWave[] = (ratioWave_log2[p] >= 1 && abs(allTwave[p] <= 0.05)) ? 3 : colorWave[p]
	colorWave[] = (ratioWave_log2[p] <= -1 && abs(allTwave[p] <= 0.05)) ? 2 : colorWave[p]
	colorWave[] = (ratioWave_log2[p] >= 1 && abs(allTwave[p] > 0.05)) ? 1 : colorWave[p]
End

STATIC Function TransformImputeBaseVal(m0,baseVal)
	WAVE m0
	Variable baseVal
	// values from Perseus - working on each replicate not on whole matrix 
	Variable width = 0.3
	Variable downShift = 1.8
	
	// make a copy of the matrix
	Duplicate/O/FREE m0,m1
	// delete base value
	if(numtype(baseVal) == 2)
		// do nothing
	else
		m1[][] = (m1[p][q] == baseVal) ? NaN : m1[p][q]
	endif
	// log tranform
	m1[][] = log(m1[p][q])
	
	Variable nCols = dimsize(m1,1)
	Variable meanVar,sdVar
	Variable i
	
	for(i = 0; i < nCols; i += 1)
		MatrixOp/O/FREE tempW = col(m1,i)
		WaveTransform zapnans tempW
		sdVar = sqrt(variance(tempW))
		meanVar = mean(tempW) - (sdVar * 1.8)
		sdVar = sdVar * width
		// add noise to base values in m0 and log transform real values col by col
		if(numtype(baseVal) == 2)
			// deal with NaN
			m0[][i] = (numtype(m0[p][i]) == 2) ? meanVar + gnoise(sdVar) : log(m0[p][i])
		else
			m0[][i] = (m0[p][i] == baseVal) ? meanVar + gnoise(sdVar) : log(m0[p][i])
		endif
	endfor
End

Function MakeVPlot()
	WAVE/Z allTWave,ratioWave_log2,colorWave,colorTableWave
	WAVE/Z/T volcanoLabelWave
	KillWindow/Z volcanoPlot
	Display/N=volcanoPlot/W=(35,45,430,734) allTWave vs ratioWave_log2
	SetAxis/W=volcanoPlot/A/R/N=1 left
	ModifyGraph/W=volcanoPlot log(left)=1
	Variable maxVar = max(wavemax(ratioWave_log2),abs(wavemin(ratioWave_log2)))
	Variable minPVar = wavemin(allTWave)
	minPVar = 10 ^ (floor((log(minPVar))))
	SetAxis/W=volcanoPlot bottom -maxVar,maxVar
	ModifyGraph/W=volcanoPlot mode=3,marker=19,mrkThick=0
	SetDrawEnv/W=volcanoPlot xcoord= bottom,ycoord= left,dash= 3;DelayUpdate
	DrawLine/W=volcanoPlot -maxVar,0.05,maxVar,0.05
	SetDrawEnv/W=volcanoPlot xcoord= bottom,ycoord= left,dash= 3;DelayUpdate
	DrawLine/W=volcanoPlot 0,1,0,minPVar
	ModifyGraph/W=VolcanoPlot zColor(allTWave)={colorWave,0,3,cindexRGB,0,colorTableWave}
	Label/W=VolcanoPlot left "P-Value"
	String labelStr = volcanoLabelWave[0] + " / " + volcanoLabelWave[1] + " (Log\\B2\\M)"
	Label/W=VolcanoPlot bottom labelStr
	SetWindow VolcanoPlot, hook(modified)=thunk_hook
End

Function TableInterestingValues()
	WAVE colorWave, allTWave, ratioWave, ratioWave_log2
	WAVE/T SHORTNAME, NAME // names of proteins hardcoded here
	Duplicate/O allTWave, allTWave_log10
	allTWave_log10 = -log(allTWave[p])
	// manhattan distance 
	MatrixOp/O productWave = abs(allTWave_log10) + abs(ratioWave_log2)
	Duplicate/O allTWave, so_allTWave
	Duplicate/O ratioWave, so_ratioWave
	Duplicate/O productWave, so_productWave
	Duplicate/O colorWave, so_colorWave
	Duplicate/O NAME, so_NAME
	Duplicate/O SHORTNAME, so_SHORTNAME
	
	Make/O/N=(numpnts(ratioWave)) keyW=p
	Duplicate/O keyW,so_keyW
	
	Sort/R {so_colorWave,so_productWave}, so_allTWave, so_ratioWave, so_productWave, so_colorWave, so_NAME, so_SHORTNAME, so_keyW
	KillWindow/Z rankTable
	Edit/N=rankTable/W=(432,45,942,734) so_NAME, so_SHORTNAME, so_productWave, so_colorWave, so_allTWave, so_ratioWave, so_keyW
End

Function AddSignificantHitsToVolcano()
	WAVE/Z so_colorWave
	WAVE/Z/T so_SHORTNAME
	String labelStr = "\Z09"
	
	Variable nRows = numpnts(so_ColorWave)
	Variable i
	
	for(i = 0; i < nRows; i += 1)
		if(so_colorWave[i] == 3)
			if(strlen(so_SHORTNAME[i]) > 0)
				labelStr += so_SHORTNAME[i] + "\r"
			endif
		elseif(so_colorWave[i] == 2)
			break
		endif
	endfor
	
	TextBox/W=volcanoPlot/C/N=topProts/F=0/A=LT/X=0.00/Y=0.00 labelStr
End

Function MakeMeanComparison()
	KillWindow/Z meanPlot
	WAVE/Z meanCond1,meanCond2,colorWave,colorTableWave
	WAVE/Z/T volcanoLabelWave
	if(!WaveExists(meanCond1) || !WaveExists(meanCond2))
		Print "Error: no mean waves."
		return -1
	endif
	
	Display/N=meanPlot/W=(944,45,1340,401) meanCond1 vs meanCond2
	Variable maxVar = max(wavemax(meanCond1),wavemax(meanCond2))
	Variable minVar = min(wavemin(meanCond1),wavemin(meanCond2))
		maxVar = 10 ^ (ceil((log(maxVar))))
		minVar = 10 ^ (floor((log(minVar)))) // if any Nans will this fail?
	SetAxis/W=meanPlot left minVar,maxVar
	SetAxis/W=meanPlot bottom minVar,maxVar
	ModifyGraph/W=meanPlot log=1
	ModifyGraph/W=meanPlot width={Plan,1,bottom,left}
	ModifyGraph/W=meanPlot mode=3,marker=19
	ModifyGraph/W=meanPlot zColor(meanCond1)={colorWave,0,3,cindexRGB,0,colorTableWave}
	ModifyGraph/W=meanPlot msize=3,mrkThick=0
	Label/W=meanPlot left volcanoLabelWave[0]
	Label/W=meanPlot bottom volcanoLabelWave[1]
	SetWindow meanPlot, hook(modified)=thunk_hook
End

STATIC Function FromVolcanoToPCA()
	WAVE/Z allCond1, allCond2
	// we will use imputed values
	Concatenate/O/NP=1 {allCond1,allCond2}, forPCA
	DoThePCA()
End

STATIC Function GetReadyForPCA(STRING selectedWavesList,VARIABLE baseVal)
	Concatenate/O selectedWavesList, forPCA
	// deal with baseVal
	TransformImputeBaseVal(forPCA,baseVal)
	DoThePCA()
End

Function DoThePCA()
	KillWindow/Z pcaPlot
	WAVE/Z forPCA
	if(!WaveExists(forPCA))
		DoAlert 0, "Missing a 2D wave, forPCA"
		return 0
	endif
	Variable nProt = dimsize(forPCA,0)
	
	WAVE/Z colorWave, colorTableWave
	if(!WaveExists(colorWave) || nProt != DimSize(forPCA,0))
		Make/O/N=(nProt) colorWave=0
		MakeColorTableWave()
		WAVE/Z colorTableWave
	endif
	if(!WaveExists(colorTableWave))
		MakeColorTableWave()
		WAVE/Z colorTableWave
	endif
	// pre-process matrix, centre on 0, sd of 1
	// we're interested in Rows (proteins) not Cols (expts)
	MatrixOp/O forPCA = SubtractMean(forPCA,2)
	MatrixOp/O forPCA = NormalizeRows(forPCA)
	// abfter subtraction we can get 0 across the row, gives NaN after norm so
	forPCA[][] = (numtype(forPCA[p][q]) == 2) ? 0 : forPCA[p][q]
	// do the PCA, SRMT flag is needed to get M_R
	PCA/ALL/SRMT forPCA
	WAVE/Z M_R
	// make a version of colorwave to display interesting proteins
	WAVE/Z allTWave
	Duplicate/O colorWave, colorPCAWave
	colorPCAWave[] = (abs(allTwave[p] > 0.05)) ? NaN : colorWave[p]
	// display PC1 and PC2
	Display/N=pcaPlot/W=(36,757,431,965) M_R[][1] vs M_R[][0]
	SetAxis/W=pcaPlot left -1,1
	SetAxis/W=pcaPlot bottom -1,1
	ModifyGraph/W=pcaPlot mode=3,marker=19,msize=2,mrkThick=0
	ModifyGraph/W=pcaPlot zColor(M_R)={colorPCAWave,0,3,ctableRGB,0,colorTableWave}
	ModifyGraph/W=pcaPlot zero=4,mirror=1
	Label/W=pcaPlot left "PC2"
	Label/W=pcaPlot bottom "PC1"
	ModifyGraph/W=pcaPlot height={Plan,1,left,bottom}
	SetWindow pcaPlot, hook(modified)=thunk_hook
End

Function LabelTopXProts(numProteins)
	Variable numProteins
	if(strlen(WinList("volcanoPlot",";","WIN:1")) < 1)
		DoAlert 0, "Make Volcano Plot first"
		return -1
	endif
	WAVE/Z so_keyW
	Make/O/N=(numProteins) topTenPos = so_KeyW[p]
	WAVE/Z/T so_SHORTNAME
	Make/O/N=(numProteins)/T topTenLabel = so_SHORTNAME[p]
	
	Variable i
	
	for(i = 0; i < numProteins; i += 1)
		Tag/A=LB/C/N=$("top"+num2str(i))/B=1/F=0/S=3/V=1/L=0/X=1/Y=1/W=volcanoPlot allTWave, topTenPos[i], topTenLabel[i]
	endfor
End

STATIC Function MakeTheLayout()
	KillWindow/Z summaryLayout
	NewLayout/N=summaryLayout
	
	AppendLayoutObject/W=summaryLayout graph meanPlot
	AppendLayoutObject/W=summaryLayout graph pcaPlot
	AppendLayoutObject/W=summaryLayout graph volcanoPlot
	LayoutPageAction/W=summaryLayout size(-1)=(595, 842), margins(-1)=(18, 18, 18, 18)
	ModifyLayout/W=summaryLayout units=0
	ModifyLayout/W=summaryLayout frame=0,trans=1
	
	ModifyLayout/W=summaryLayout left(meanPlot)=322,top(meanPlot)=21,width(meanPlot)=258,height(meanPlot)=222
	ModifyLayout/W=summaryLayout left(pcaPlot)=322,top(pcaPlot)=244,width(pcaPlot)=260,height(pcaPlot)=224
	ModifyLayout/W=summaryLayout left(volcanoPlot)=21,top(volcanoPlot)=21,height(volcanoPlot)=450,width(volcanoPlot)=320
End


////////////////////////////////////////////////////////////////////////
// Panel functions
////////////////////////////////////////////////////////////////////////

Function VolcanoIO_Panel()
	WAVE/Z/T volcanoPrefixWave, volcanoLabelWave
	if(!WaveExists(volcanoPrefixWave))
		MAKE/O/N=(2)/T volcanoPrefixWave = {"prefix1*","prefix2*"}
		MAKE/O/N=(2)/T volcanoLabelWave = {"Test","Control"}
	endif
	String prefix1 = volcanoPrefixWave[0]
	String prefix2 = volcanoPrefixWave[1]
	String label1 = volcanoLabelWave[0]
	String label2 = volcanoLabelWave[1]
	Variable baseVal = 0
	Variable pairOpt = 0
	
	DoWindow/K VolcanoSetup
	NewPanel/N=VolcanoSetup/K=1/W=(81,73,774,200)
	SetVariable box1,pos={76,13},size={200,14},title="Prefix for condition 1 (test):",value=_STR:prefix1
	SetVariable box2,pos={76,44},size={200,14},title="Prefix for condition 2 (ctrl):",value=_STR:prefix2
	SetVariable box3,pos={276,13},size={200,14},title="Label for condition 1:",value=_STR:label1
	SetVariable box4,pos={276,44},size={200,14},title="Label for condition 2:",value=_STR:label2
	SetVariable box5,pos={19,75},size={250,14},title="What value represents absent proteins?",format="%g",value=_NUM:baseVal
	CheckBox box6,pos={208,106},size={20,20},title="Analyse paired data?",value=pairOpt,mode=0
	Button DoIt,pos={564,100},size={100,20},proc=ButtonProc,title="Do It"
End
 
Function ButtonProc(ba) 
	STRUCT WMButtonAction &ba
	if( ba.eventCode != 2 )
		return 0
	endif
	ControlInfo/W=$ba.win box1
	String prefix1 = S_Value
	ControlInfo/W=$ba.win box2
	String prefix2 = S_Value
	ControlInfo/W=$ba.win box3
	String label1 = S_Value
	ControlInfo/W=$ba.win box4
	String label2 = S_Value
	ControlInfo/W=$ba.win box5
	Variable baseVal = V_Value
	ControlInfo/W=$ba.win box6
	Variable pairOpt = V_Value
	Print "Test group:", prefix1, "Labelled", label1, "\rControl group:", prefix2, "Labelled", label1, "\rValue for imputation:", baseVal
	if(pairOpt == 1)
		Print "Using paired data."
	endif
	WAVE/Z/T 	volcanoPrefixWave, volcanoLabelWave
	volcanoPrefixWave[0] = prefix1
	volcanoPrefixWave[1] = prefix2
	volcanoLabelWave[0] = label1
	volcanoLabelWave[1] = label2
	VolcanoWorkflowWrapper(prefix1,prefix2,baseVal,pairOpt)
End

Function MakePCAWaveSelectorPanel()
	
	String panelName = "PCASelector"
	Variable baseVal = 0
	
	if (WinType(panelName) == 7)
		// if the panel already exists, show it
		DoWindow/F $panelName
	else
		// doesn't exist, make it
		NewPanel/K=1/N=$panelName/W=(181,179,471,540) as "Select Waves for PCA"
		// list box control doesn't have any attributes set on it
		ListBox ExampleWaveSelectorList,pos={9,13},size={273,241}
		// This function does all the work of making the listbox control into a
		// Wave Selector widget. Note the optional parameter that says what type of objects to
		// display in the list. 
		MakeListIntoWaveSelector(panelName, "ExampleWaveSelectorList", content = WMWS_Waves)

		// This function does all the work of making a PopupMenu control into a wave sorting control
		PopupMenu sortKind, pos={9,270},title="Sort Waves By"
		MakePopupIntoWaveSelectorSort(panelName, "ExampleWaveSelectorList", "sortKind")
		SetVariable box3,pos={9,300},size={230,14},title="What value represents absent proteins?",format="%g",value=_NUM:baseVal
		Button doPCA,pos={9,330},size={110,20},proc=doPCAButtonProc,title="Do PCA"
	endif
End

Function doPCAButtonProc(ctrlName) : ButtonControl
	String ctrlName
	
	ControlInfo/W=PCASelector box3
	Variable baseVal = V_Value
	String selectedWavesList = WS_SelectedObjectsList("PCASelector","ExampleWaveSelectorList")
	GetReadyForPCA(selectedWavesList,baseVal)
End

////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

// Modified from _sk http://www.igorexchange.com/node/7797
Function thunk_hook(s)
	Struct WMWinHookStruct& s
	WAVE allTWave
	WAVE/T SHORTNAME
 
	strswitch (s.eventname)
		case "mouseup":

			String s_traceinfo = TraceFromPixel(s.mouseloc.h, s.mouseloc.v, "WINDOW:"+s.winname+";")
			Variable v_pt = str2num(StringByKey("HITPOINT", s_traceinfo))
			String targetTrace = StringByKey("TRACE", s_traceinfo) // now takes whatever trace is clicked on
 
			if (numtype(v_pt) != 2)
				Tag/c/n=t1/b=3/f=0/s=3/v=1/X=10/Y=10 $targetTrace, v_pt, SHORTNAME[v_pt]
			else
				Tag/n=t1/k
			endif
			break
 
		case "kill":
			setwindow $(s.winname), hook(modified)=$""
			break
	endswitch
 
end

STATIC Function MakeColorTableWave()
	Make/O/N=(4,4) colorTableWave = 32768
	colorTableWave[1][0] = 65535
	colorTableWave[3][0] = 65535
	colorTableWave[2][2] = 65535
	colorTableWave[3][2] = 65535
End