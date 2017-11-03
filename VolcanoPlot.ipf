#pragma TextEncoding = "MacRoman"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Math Utility Functions>

Menu "Macros"
	"Volcano Plot...",  VolcanoIO_Panel()
End

Function VolcanoIO_Panel()
	
	String prefix1 = "prefix1*"
	String prefix2 = "prefix2*"
	Variable baseVal = 1
	
	DoWindow/K VolcanoSetup
	NewPanel/N=VolcanoSetup/K=1/W=(81,73,774,200)
	SetVariable box1,pos={86,13},size={500,14},title="Prefix for condition 1",value=_STR:prefix1
	SetVariable box2,pos={86,44},size={500,14},title="Prefix for condition 2",value=_STR:prefix2
	SetVariable box3,pos={86,75},size={166,14},title="Value for absent proteins",format="%g",value=_NUM:baseVal
	
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
	Concatenate/O wList1, allCond1
	String wList2 = WaveList(prefix2,";","")
	Concatenate/O wList2, allCond2
	
	// deal with baseVal
	ImputeBaseVal(allCond1,baseVal)
	ImputeBaseVal(allCond2,baseVal)
	
	// make mean waves
	MatrixOp/O meanCond1 = sumrows(allCond1)
	meanCond1 /=ItemsInList(wList1)
	MatrixOp/O meanCond2 = sumrows(allCond2)
	meanCond2 /=ItemsInList(wList2)
	
	// ratio wave
	MatrixOp/O ratioWave = meanCond1 / meanCond2

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
	Duplicate/O ratioWave,ratioWave_log2
	ratioWave_log2 = log2(abs(ratioWave[p]))
	// assign colors
	colorWave = (abs(ratioWave_log2 >= 1)) ? colorWave[p] + 1 : colorWave[p]
	colorWave = (abs(allTwave < 0.05)) ? colorWave[p] + 2 : colorWave[p]
	MakeColorTableWave()
	MakeVPlot()
	TableInterestingValues()
End

STATIC Function ImputeBaseVal(m0,baseVal)
	WAVE m0
	Variable baseVal
	
	// make a copy of the matrix
	Duplicate/O/FREE m0,m1
	// delete base value and then find the mean and sd of lowest 15 values
	if(numtype(baseVal) == 2)
		// do nothing
	else
		m1[][] = (m1[p][q] == baseVal) ? NaN : m1[p][q]
	endif
	Redimension/N=(dimsize(m0,0)*dimsize(m0,1)) m1
	WaveTransform zapnans m1
	Sort m1, m1
	
	Variable meanVar = mean(m1,0,14)
	Variable sdVar = sqrt( variance(m1,0,14) )
	
	// add noise to base values in m0
	if(numtype(baseVal) == 2)
		// deal with NaN
		m0[][] = (numtype(m0[p][q]) == 2) ? meanVar + gnoise(sdVar) : m0[p][q]
	else
		m0[][] = (m0[p][q] == baseVal) ? meanVar + gnoise(sdVar) : m0[p][q]
	endif
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
	ModifyGraph/W=volcanoPlot mode=3,marker=19
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

// from _sk http://www.igorexchange.com/node/7797

function thunk_make_data()
	string s_winname = "thunkdata"
 
	dowindow $s_winname
	if (v_flag)
		dowindow/k $s_winname
	endif
 
	make/o/n=10 w_a = gnoise(10)
	make/o/n=10 w_b = gnoise(10)
	make/o/n=10 w_c = gnoise(10)
	display/n=$s_winname/k=1 w_a vs w_b
	ModifyGraph/w=$s_winname mode=3,marker=16,msize=5,useMrkStrokeRGB=1;DelayUpdate
	ModifyGraph/w=$s_winname zColor(w_a)={w_c,*,*,Grays,0};DelayUpdate
	ModifyGraph/w=$s_winname width={Aspect,1},height={Aspect,1};DelayUpdate
	ModifyGraph/w=$s_winname mirror=2
 
	setwindow $s_winname, hook(modified)=thunk_hook
end
 
function thunk_hook(s)
	struct WMWinHookStruct& s
	WAVE allTWave
	WAVE/T SHORTNAME
 
	strswitch (s.eventname)
		case "mouseup":
			wave w_c = root:w_c
			wave w_a = root:w_a
			wave w_b = root:w_b
 
			string s_traceinfo = TraceFromPixel(s.mouseloc.h, s.mouseloc.v, "WINDOW:"+s.winname+";")
			variable v_pt = str2num(stringbykey("HITPOINT", s_traceinfo))
 
			if (numtype(v_pt) != 2)
				//print "f(z): ", w_c[v_pt]
				tag/c/n=t1/b=3/f=0/s=3/v=1/X=10/Y=10 allTwave, v_pt, SHORTNAME[v_pt]
			else
				tag/n=t1/k
			endif
			break
 
		case "kill":
			setwindow $(s.winname), hook(modified)=$""
			break
	endswitch
 
end

