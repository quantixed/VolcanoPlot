#pragma TextEncoding = "MacRoman"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Math Utility Functions>

// Figure in the paper was generated with
// https://github.com/quantixed/TPD54/commit/accb6a86619aa4668e5dadac655e55630e9a55f2
// Load data into Igor and then run the command from the Macro menu.
// Note that naming of the columns is important - check example file for details.

Menu "Macros"
	"Volcano Plot...",  VolcanoIO_Panel()
End

Function VolcanoIO_Panel()
	
	String prefix1 = "prefix1*"
	String prefix2 = "prefix2*"
	Variable baseVal = 0
	
	DoWindow/K VolcanoSetup
	NewPanel/N=VolcanoSetup/K=1/W=(81,73,774,200)
	SetVariable box1,pos={86,13},size={500,14},title="Prefix for condition 1 (test)",value=_STR:prefix1
	SetVariable box2,pos={86,44},size={500,14},title="Prefix for condition 2 (ctrl)",value=_STR:prefix2
	SetVariable box3,pos={29,75},size={300,14},title="What value represents absent proteins?",format="%g",value=_NUM:baseVal
	
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
    Print "Test group is", prefix1, "Control group is", prefix2, "Value for imputation is", V_Value
    MakeVolcano(prefix1,prefix2,V_Value)
End

// This function drives the whole program
/// @param	prefix1	string prefix that will select all waves of condition1
/// @param	prefix2	string prefix that will select all waves of condition1
/// @param	baseVal	proteins that are absent have this value
Function MakeVolcano(prefix1,prefix2,baseVal)
	String prefix1,prefix2
	Variable baseVal
	
	String wList1 = WaveList(prefix1,";","")
	String wList2 = WaveList(prefix2,";","")
	if(ItemsInList(wList1) == 0 || ItemsInList(wList2) == 0)
		DoAlert 0, "Missing data"
		return -1
	endif
	Concatenate/O wList1, allCond1
	Concatenate/O wList2, allCond2
	
	// deal with baseVal
	TransformImputeBaseVal(allCond1,baseVal)
	TransformImputeBaseVal(allCond2,baseVal)
	
	// now do T-tests
	Variable nProt = dimsize(allCond1,0)
	Make/O/N=(nProt) allTWave,colorWave=0
	Variable pVar
	WAVE/Z W_StatsTTest
	
	Variable i
	
	for(i = 0; i < nProt; i += 1)
		MatrixOp/O/FREE w0 = row(allCond1,i) ^ t
		MatrixOp/O/FREE w1 = row(allCond2,i) ^ t
		StatsTTest/Q/Z w0,w1
		WAVE/Z W_StatsTTest
		if(V_flag == 0)
			pVar = W_StatsTTest[9] // p-value
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
	ratioWave_log2 = log2(abs(ratioWave[p]))
	// assign colors
	colorWave = (abs(ratioWave_log2 >= 1)) ? colorWave[p] + 1 : colorWave[p]
	colorWave = (abs(allTwave < 0.05)) ? colorWave[p] + 2 : colorWave[p]
	MakeColorTableWave()
	MakeVPlot()
	TableInterestingValues()
	MakeMeanComparison()
	DoThePCA()
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
	WAVE allTWave,ratioWave_log2,colorWave,colorTableWave
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
	SetWindow VolcanoPlot, hook(modified)=thunk_hook
End

STATIC Function MakeColorTableWave()
	Make/O/N=(4,4) colorTableWave = 32768
	colorTableWave[1][0] = 65535
	colorTableWave[3][0] = 65535
	colorTableWave[2][2] = 65535
	colorTableWave[3][2] = 65535
End

Function TableInterestingValues()
	WAVE colorWave,allTWave,ratioWave
	WAVE/T SHORTNAME,NAME // names of proteins hardcoded here - maybe add this to the panel?
	Duplicate/O allTWave, allTWave_log10
	allTWave_log10 = -log(allTWave[p])
	MatrixOp/O productWave = allTWave_log10 * ratioWave
	Duplicate/O allTWave, so_allTWave
	Duplicate/O ratioWave, so_ratioWave
	Duplicate/O productWave, so_productWave
	so_productWave = abs(productWave[p])
	Duplicate/O colorWave, so_colorWave
	Duplicate/O NAME, so_NAME
	Duplicate/O SHORTNAME, so_SHORTNAME
	
	Make/O/N=(numpnts(ratioWave)) keyW=p
	Duplicate/O keyW,so_keyW
	
	Sort/R {so_colorWave,so_productWave}, so_allTWave, so_ratioWave, so_productWave, so_colorWave, so_NAME, so_SHORTNAME, so_keyW
	KillWindow/Z rankTable
	Edit/N=rankTable/W=(432,45,942,734) so_NAME, so_SHORTNAME, so_productWave, so_colorWave, so_allTWave, so_ratioWave, so_keyW
End

Function MakeMeanComparison()
	KillWindow/Z meanPlot
	WAVE/Z meanCond1,meanCond2
	WAVE/Z colorWave,colorTableWave
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
	Label/W=meanPlot left "Condition 1"
	Label/W=meanPlot bottom "Condition 2"
	SetWindow meanPlot, hook(modified)=thunk_hook
End

Function DoThePCA()
	KillWindow/Z pcaPlot
	WAVE/Z colorWave, colorTableWave, allCond1, allCond2
	Concatenate/O {allCond1,allCond2}, forPCA
	PCA/ALL/SRMT forPCA
	WAVE/Z M_R
	Display/N=pcaPlot/W=(36,757,431,965) M_R[][1] vs M_R[][0]
	WaveStats/RMD=[][0,1]/Q M_R
	SetAxis/W=pcaPlot left V_min,V_Max
	SetAxis/W=pcaPlot bottom V_min,V_Max
	ModifyGraph/W=pcaPlot mode=3,marker=19,mrkThick=0,zColor(M_R)={colorWave,*,*,ctableRGB,0,colorTableWave}
	ModifyGraph/W=pcaPlot zero=4,mirror=1
	Label/W=pcaPlot left "PC2"
	Label/W=pcaPlot bottom "PC1"
	ModifyGraph/W=pcaPlot width={Plan,1,bottom,left}
	SetWindow meanPlot, hook(modified)=thunk_hook
End

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