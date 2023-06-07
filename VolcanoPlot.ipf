#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later
#include <WaveSelectorWidget>
#include <PopupWaveSelector>

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////
Menu "Macros"
	SubMenu "Proteomics"
		"Load MaxQuant Data...", /Q, LoadMaxQuantData()
		"Volcano Plot...", /Q, VolcanoIO_Panel()
		SubMenu "Subcellular analysis"
			"Make List to Retrieve Uniprot Data", /Q, UniprotTable()
			"Load and Match UniProt Data...", /Q, UniprotWrapper()
			"Filter for GO Term(s)", /Q, GOTerm_Panel()
		end
//		"PCA Only...", /Q, MakePCAWaveSelectorPanel() // remove this option
		"Label Top 10", /Q, LabelTopTenWorkflow()
		"Count Comparison", /Q, MakeTheComparison()
		"Save Layout", /Q, SaveTheLayout()
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

Function VolcanoWorkflowWrapper()
	MakeVolcano()
	MakeColorTableWave()
	MakeVPlot()
	TableInterestingValues()
	AddSignificantHitsToVolcano(0)
	MakeMeanComparison()
	FromVolcanoToPCA()
	MakeTheLayout()
//	SaveTheLayout()
End

Function UniprotWrapper()
	if(LoadUniprot() == 0)
		MatchResultsAndUniProt()
	endif
	if(FilterGOTerms() == 0)
		// build listbox that allows multiple selections
		GOTerm_Panel()
		// this panel triggers FilteredVolcanoWorkflowWrapper()
	endif
	// use selections to filter gene-names that feature the GO terms
	// plot volcano using only those gene_names
End

Function FilteredVolcanoWorkflowWrapper()
	FilterForGOTerms()
	MakeFilteredVPlot()
	AddSignificantHitsToVolcano(1)
	AddGOTermsToVolcano()
End

Function LabelTopTenWorkflow()
	LabelTopXProts(10)
End

////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////

Function LoadMaxQuantFile()
	NewDataFolder/O/S root:data 
	LoadWave/A/D/J/K=0/L={0,0,0,0,0}/V={"\t","",0,0}/W/O/Q ""
	if (strlen(S_Path) == 0) // user pressed cancel
		return -1
	endif
	NewPath/O/Q DiskFolder, S_Path
	WAVE/Z/T Gene_names, Protein_names
	Duplicate/O Gene_names, root:SHORTNAME
	Duplicate/O Protein_names, root:NAME
	Wave/T SHORTNAME = $("root:SHORTNAME")
	Wave/T NAME = $("root:NAME")
	String wList = WaveList("LFQ_Intensity*",";","")
	Variable nWaves = ItemsInList(wList)
	String newList = "root:SHORTNAME;root:NAME;"
	String wName, newName
	
	Variable i,j
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		newName = "root:" + ReplaceString("LFQ_Intensity_",wName,"")
		Duplicate/O $wName, $newName
		newList += newName + ";"
	endfor
	
	// remove data where there is a blank in both SHORTNAME and NAME
	nWaves = ItemsInList(newList)
	Variable nRow = numpnts(Gene_Names)
	
	for(i = nRow - 1; i >= 0; i -= 1)
		if(strlen(Gene_Names[i]) == 0 && strlen(Protein_names[i]) == 0)
			for(j = 0; j < nWaves; j += 1)
				DeletePoints i,1,$(StringFromList(j,newList))
			endfor
		endif
	endfor
	
	SetDataFolder root:
	return 0
End

// This function drives the whole program
// Params are stored in two waves in root volcanoPrefixWave and volcanoParamWave
Function MakeVolcano()
	WAVE/Z volcanoParamWave
	WAVE/Z/T volcanoPrefixWave
	String prefix1 = volcanoPrefixWave[0] // string prefix that will select all waves of condition 1
	String prefix2 = volcanoPrefixWave[1] // string prefix that will select all waves of condition 2
	Variable baseVal = volcanoParamWave[0] // proteins that are absent have this value
	Variable pairOpt = volcanoParamWave[1] // 1 is paired, 0 is not
	Variable seedOpt = volcanoParamWave[2] // 1 is use seed, 0 is not (truly random)
	Variable foldChange = log(abs(volcanoParamWave[3])) / log(2) // log2 value for threshold i.e. 1 is 2-fold change
	Variable seed1 = volcanoParamWave[4] // 1-1000 value for seed
	Variable seed2 = volcanoParamWave[5] // 1-1000 value for seed
	
	String wList1 = WaveList(prefix1,";","")
	String wList2 = WaveList(prefix2,";","")
	if(ItemsInList(wList1) == 0 || ItemsInList(wList2) == 0)
		DoAlert 0, "Missing data"
		return -1
	else
		Print "Test runs:", ItemsInList(wList1), "\rControl runs:", ItemsInList(wList2)
	endif
	// make sure the lists are in order
	wList1 = SortList(wList1)
	wList2 = SortList(wList2)
	// Possibly we should use a more sophisticated way to check the groups match?
	// using this method x_1, x_2, x_3 will run with y_1, y_4, y_5
	Concatenate/O/NP=1 wList1, allCond1
	Concatenate/O/NP=1 wList2, allCond2
	
	// deal with baseVal
	TransformImputeBaseVal(allCond1)
	TransformImputeBaseVal(allCond2)
	
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
				pVar = W_StatsTTest[%P] // p-value
			else
				pVar = W_StatsTTest[%P] // p-value for Paired
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
	ratioWave_log2[] = log(abs(ratioWave[p])) / log(2)
	// assign colors
	colorWave[] = (ratioWave_log2[p] >= foldChange && abs(allTwave[p] <= 0.05)) ? 3 : colorWave[p]
	colorWave[] = (ratioWave_log2[p] <= -foldChange && abs(allTwave[p] <= 0.05)) ? 2 : colorWave[p]
	colorWave[] = (ratioWave_log2[p] >= foldChange && abs(allTwave[p] > 0.05)) ? 1 : colorWave[p]
End

/// @param m0 2D wave containing data to be imputed
STATIC Function TransformImputeBaseVal(m0)
	WAVE m0
	WAVE/Z volcanoParamWave
	// work out which wave we are doing condition 1 or condition 2
	String wName = NameOfWave(m0)
	Variable nn = str2num(wName[strlen(wName) - 1])
	Variable baseVal = volcanoParamWave[0] 
	Variable seed = volcanoParamWave[3 + nn] / 1000 // will be row 4 for condition 1 and 5 for 2
	// values from Perseus - working on each replicate not on whole matrix 
	Variable width = 0.3
	Variable downShift = 1.8
	if(volcanoParamWave[2] != 0)
		SetRandomSeed seed
	endif
	
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

Function AddSignificantHitsToVolcano(vpOpt)
	Variable vpOpt // 0 is original VP, 1 is filtered VP
	
	String plotName
	if(vpOpt == 0)
		Wave/Z colorW = so_colorWave
		Wave/Z/T shortnameW = so_SHORTNAME
		plotName = "volcanoPlot"
	else
		Wave/Z colorW = ft_colorWave
		Wave/Z/T shortnameW = ft_SHORTNAME
		plotName = "ftVolcanoPlot"
	endif
	
	String labelStr = "\Z09"
	
	Variable nRows = numpnts(colorW)
	Variable i
	
	for(i = 0; i < nRows; i += 1)
		if(colorW[i] == 3)
			if(strlen(shortnameW[i]) > 0)
				labelStr += shortnameW[i] + "\r"
			endif
		elseif(colorW[i] == 2)
			break
		endif
	endfor
	
	TextBox/W=$plotName/C/N=topProts/B=1/F=0/A=LT/X=0.00/Y=0.00 labelStr
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
	WAVE/Z forPCA
	KillWaves/Z forPCA
	// we will use imputed values and we will take the ratio per replicate
	// note that Volcano plot ratio is the mean of cond1 / mean of cond2
	// here we take the ratio of replicates and use the replicates for PCA, so...
	if(DimSize(allCond1,1) != DimSize(allCond2,1))
		Print "Number of replicates between conditions differs. PCA is of all data, not ratios."
		Concatenate/O/NP=1 {allCond1,allCond2}, forPCA
	elseif(DimSize(allCond1,1) == 1)
		Print "Single replicate dataset. No PCA."
		Concatenate/O/NP=1 {allCond1,allCond2}, forPCA
		return 0
	else
		MatrixOp/O forPCA = allCond1 / allCond2
	endif
	DoThePCA()
End

STATIC Function GetReadyForPCA(STRING selectedWavesList,VARIABLE baseVal)
	Concatenate/O selectedWavesList, forPCA
	// deal with baseVal - this uses random values even if a seed was used for other analyses
	TransformImputeBaseVal(forPCA)
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
	// pre-process data
	// Previously we did
//	MatrixOp/O forPCA = SubtractMean(forPCA,2)
//	MatrixOp/O forPCA = NormalizeRows(forPCA)
	// now we do:
	// log2 transform
	MatrixOp/O forPCA = log2(forPCA)
	// in case we have NaN or Inf
	forPCA[][] = (numtype(forPCA[p][q]) == 2 || numtype(forPCA[p][q]) == 1) ? 0 : forPCA[p][q]
	// centre the data
	MatrixOp/O forPCA = SubtractMean(forPCA,2)
	// do the PCA, SRMT flag is needed to get M_R
	PCA/ALL/SRMT forPCA
	WAVE/Z M_R
	// M_R is inverse compared to the output from SIMCA-P+
	M_R *= -1
	// make a version of colorwave to display interesting proteins
	Duplicate/O colorWave, colorPCAWave
	// display PC1 and PC2
	Display/N=pcaPlot/W=(36,757,431,965) M_R[][1] vs M_R[][0]
	// find min and max for PC1 and PC2 combined
	WaveStats/Q/RMD=[][0,1] forPCA
	SetAxis/W=pcaPlot left V_min,V_max
	SetAxis/W=pcaPlot bottom V_min,V_max
	ModifyGraph/W=pcaPlot mode=3,marker=19,mrkThick=0
	ModifyGraph/W=pcaPlot zColor(M_R)={colorPCAWave,0,3,ctableRGB,0,colorTableWave}
	ModifyGraph/W=pcaPlot zmrkSize(M_R)={colorPCAWave,0,1,1,2}
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

Function MakeTheComparison()
	WAVE/Z/T volcanoPrefixWave, volcanoLabelWave
	WAVE/Z colorWave, colorTableWave
	Duplicate/O colorWave, highlightWave
	// use original values (no imputation)
	Concatenate/O/NP=1 WaveList(volcanoPrefixWave[0],";",""), countCond1
	Concatenate/O/NP=1 WaveList(volcanoPrefixWave[1],";",""), countCond2
	Variable nCol1 = DimSize(countCond1,1)
	Variable nCol2 = DimSize(countCond2,1)
	MatrixOp/O countCond1 = sumRows(countCond1) / nCol1
	MatrixOp/O countCond2 = sumRows(countCond2) / nCol2
	Variable maxVal = max(WaveMax(countCond1),WaveMax(countCond2))
	KillWindow/Z compPlot
	Display/N=compPlot/W=(436,757,700,965) countCond1 vs countCond2
	SetAxis/W=compPlot left 0,maxVal
	SetAxis/W=compPlot bottom 0,maxVal
	ModifyGraph/W=compPlot mode=3,marker=19,msize=2,mrkThick=0
	Label/W=compPlot left VolcanoLabelWave[0]
	Label/W=compPlot bottom VolcanoLabelWave[1]
	ModifyGraph/W=compPlot width={Aspect,1}
	ModifyGraph/W=compPlot zColor(countCond1)={highlightWave,0,3,ctableRGB,0,colorTableWave}
	SetWindow compPlot, hook(modified)=thunk_hook
End

// use a semi-colon separated list of protein names to highlight
Function HighlightTheseProteins(myList)
	String myList
	
	WAVE/Z highlightWave
	if(!WaveExists(highlightWave))
		WAVE/Z colorWave
		Duplicate/O colorWave, highlightWave
	endif
	
	WAVE/Z/T SHORTNAME
	Variable nProteins = ItemsInList(myList)
	highlightWave = 0
	String proteinName
	
	Variable i
	
	for(i = 0; i < nProteins; i += 1)
		proteinName = StringFromList(i, myList)
		FindValue/TEXT=proteinName/TXOP=2 SHORTNAME
		highlightWave[V_Value] = 3
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

Function SaveTheLayout()
	WAVE/Z/T volcanoLabelWave
	String fileName = "summaryLayout_" + volcanoLabelWave[0] + "vs" + volcanoLabelWave[1] + ".pdf"
	SavePICT/O/WIN=summaryLayout/P=DiskFolder/E=-2/W=(0,0,0,0) as fileName
	fileName = "summaryLayout_" + volcanoLabelWave[0] + "vs" + volcanoLabelWave[1] + ".png"
	SavePICT/O/WIN=summaryLayout/P=DiskFolder/E=-5/RES=300 as fileName
End

// This function will load an output from UniProt.
// Using a list of gene names from VolcanoPlot output (via copy/paste), query in UniProt ID mapping
// select columns for Gene ontology (cellular component) & Subcellular location [CC]
Function LoadUniprot()
	if(!DataFolderExists("root:data"))
		NewDataFolder/O root:data
	endif
	NewDataFolder/O/S root:data:uniprot
	// load the tab-separated output downloaded from UniProt	
	LoadWave/A/D/J/K=0/L={0,0,0,0,0}/V={"\t","",0,0}/W/O/Q ""
	if (strlen(S_Path) == 0) // user pressed cancel
		return -1
	endif
	
	SetDataFolder root:
	return 0
End

Function UniprotTable()
	WAVE/Z/T SHORTNAME
	Duplicate/O/T SHORTNAME, clickMe
	DeleteEmptyCellsFromTextWave(clickMe)
	KillWindow/Z closeTableWhenFinished
	Edit/N=closeTableWhenFinished/K=1 clickMe as "Close this table when finished"
	DoAlert 0, "Click `clickMe` to select the column.\rCopy data to clipboard and use at:\rhttps://www.uniprot.org/uploadlists"
End

Function MatchResultsAndUniprot()
	// first we need to check that we have results and uniprot data
	SetDataFolder root:
	WAVE/Z/T SHORTNAME
	if(!WaveExists(SHORTNAME))
		Print "Proteomics data not loaded"
		return -1
	endif
	Wave/Z/T Gene_names = root:data:uniprot:Gene_names
	Wave/Z/T GOw = root:data:uniprot:Gene_ontology__cellular_component_
	Wave/Z/T SCw = root:data:uniprot:Subcellular_location__CC_
	if(!WaveExists(Gene_names) || !WaveExists(GOw) || !WaveExists(SCw))
		Print "Uniprot data not loaded"
		return -1
	endif
	
	Variable nHits = numpnts(SHORTNAME)
	Variable nGenes = numpnts(Gene_names)
	
	// first make a textwave of Gene_names where the alternatives are semi-colon separated
	Duplicate/O/FREE/T Gene_Names, gnW
	gnW[] = ReplaceString(" ",gnW[p],";") + ";"
	// make the wave to hold the GO terms and subcellular info
	Make/O/T/N=(nHits) GO_Terms, SC_Terms
	String sName, gNameList
	Variable matchVar
	
	Variable i,j,k
	
	for(i = 0; i < nHits; i += 1)
		sName = SHORTNAME[i]
		if(strlen(sName) == 0)
			continue
		endif
		matchVar = 0
		
		for(j = 0; j < nGenes; j += 1)
			gNameList = Gene_names[j]
			
			for(k = 0; k < ItemsInList(gNameList); k += 1)
				if(cmpstr(sName,StringFromList(k,gNameList)) == 0)
					GO_Terms[i] = GOw[j]
					SC_Terms[i] = SCw[j]
					matchVar = 1
					break
				endif
			endfor
			// if we found a match we will break and go back to i
			if(matchVar == 1)
				break
			endif
		endfor
	endfor
	
	// now map across to sorted wave if it is there. If it isn't Volcano plot won't map them as the code stands now
	WAVE/Z so_KeyW
	if(!WaveExists(so_KeyW))
		return -1
	endif
	Make/O/T/N=(numpnts(so_KeyW)) so_GO_Terms = GO_Terms[so_keyW[p]]
	Make/O/T/N=(numpnts(so_KeyW)) so_SC_Terms = SC_Terms[so_keyW[p]]
End

Function FilterGOTerms()
	WAVE/Z GO_Terms
	Duplicate/O/FREE/T GO_Terms, tempGT
	// we might have a space at the beginning of each item
//	tempGT[] = SelectString(cmpstr(" ",tempGT[p][0][0][0][0]),tempGT[p][0][0][0][1,strlen(tempGT[p])-1],tempGT[p])
	DeleteEmptyCellsFromTextWave(tempGT)
	tempGT[] = ReplaceString("; ",tempGT[p],";")
	// make list of all GO Terms. Rows have multiple GO terms separated by semi-colon
	// row does not terminate with semicolon
	String goList
	wfprintf goList, "%s;", tempGT
	// Case-insensitive alphanumeric sort
	goList = SortList(goList,";",8)
	// convert to text wave
	Wave/T bigGOw = ListToTextWave(goList,";")
	FindDuplicates/RT=all_GO_Terms bigGOw
	WAVE/Z/T all_GO_terms
	Make/O/N=(numpnts(all_GO_Terms))/B all_GO_Terms_sel = 0
	
	return 0
End

Function FilterForGOTerms()
	// first make a textwave of the GO terms we will filter for
	WAVE/Z/T all_GO_Terms
	WAVE/Z all_GO_Terms_sel
	Duplicate/O/T all_GO_Terms, filter_GO_Terms
	filter_GO_Terms[] = SelectString(all_GO_Terms_sel[p],"",all_GO_Terms[p])
	DeleteEmptyCellsFromTextWave(filter_GO_Terms)
	Variable nGOT = numpnts(filter_GO_Terms)
	
	WAVE/Z/T GO_Terms
	Variable nGenes = numpnts(GO_Terms)
	Make/O/N=(nGenes)/B GO_Terms_nFilter = 0 // zero is delete, 1 is keep
	String gostring
	
	Variable i,j
	
	for(i = 0; i < nGenes; i += 1)
		gostring = GO_Terms[i]
		
		for(j = 0; j < nGOT; j += 1)
			if(strsearch(gostring,filter_GO_Terms[j],0,2) >= 0)
				// we have a hit
				GO_Terms_nFilter[i] = 1
				break
			endif
		endfor
	endfor
	// now let's filter - we'll do it by duplicate and delete
	// we need all of these except productWave
	WAVE/Z so_allTWave, so_colorWave, so_keyW, so_productWave, so_ratioWave
	WAVE/Z/T so_GO_Terms, so_NAME, so_SC_Terms, so_SHORTNAME
	nGenes = numpnts(so_keyW)
	Duplicate/O so_allTWave, ft_allTWave
	Duplicate/O so_colorWave, ft_colorWave
	Duplicate/O so_keyW, ft_keyW
	Duplicate/O so_ratioWave, ft_ratioWave
	Duplicate/O/T so_GO_Terms, ft_GO_Terms
	Duplicate/O/T so_NAME, ft_NAME
	Duplicate/O/T so_SC_Terms, ft_SC_Terms
	Duplicate/O/T so_SHORTNAME, ft_SHORTNAME
	
	for(i = 0; i < nGenes; i += 1)
		if(GO_Terms_nFilter[so_keyW[i]] != 1) // delete
			ft_allTWave[i] = NaN
			ft_colorWave[i] = NaN
			ft_keyW[i] = NaN
			ft_ratioWave[i] = NaN
			ft_GO_Terms[i] = ""
			ft_NAME[i] = ""
			ft_SC_Terms[i] = ""
			ft_SHORTNAME[i] = ""
		else
			// duplicate means that data should be in place, however we may have blanks in Text waves
			// GO_Terms should not need checking
			if(strlen(ft_NAME[i]) == 0)
				ft_NAME[i] = "noName"
			endif
			if(strlen(ft_SC_Terms[i]) == 0)
				ft_SC_Terms[i] = "noInfo"
			endif
			if(strlen(ft_SHORTNAME[i]) == 0)
				ft_SHORTNAME[i] = "noName"
			endif
		endif
	endfor
	WaveTransform ZapNans ft_allTWave
	WaveTransform ZapNans ft_colorWave
	WaveTransform ZapNans ft_keyW
	WaveTransform ZapNans ft_ratioWave
	// need ratiowave log2 for plotting
	Duplicate/O ft_ratioWave, ft_ratioWave_log2
	ft_ratioWave_log2[] = log(abs(ft_ratioWave[p])) / log(2)
	DeleteEmptyCellsFromTextWave(ft_GO_Terms)
	DeleteEmptyCellsFromTextWave(ft_NAME)
	DeleteEmptyCellsFromTextWave(ft_SC_Terms)
	DeleteEmptyCellsFromTextWave(ft_SHORTNAME)
	KillWindow/Z ftRankTable
	Edit/N=ftRankTable/W=(432,45,942,734) ft_NAME, ft_SHORTNAME, ft_GO_Terms, FT_SC_Terms, ft_colorWave, ft_allTWave, ft_ratioWave, ft_keyW
End

Function MakeFilteredVPlot()
	String plotName = "ftVolcanoPlot"
	WAVE/Z ft_allTWave,ft_ratioWave_log2,ft_colorWave,colorTableWave
	WAVE/Z/T volcanoLabelWave
	KillWindow/Z $plotName
	Display/N=$plotName/W=(35,45,430,734) ft_allTWave vs ft_ratioWave_log2
	SetAxis/W=$plotName/A/R/N=1 left
	ModifyGraph/W=$plotName log(left)=1
	Variable maxVar = max(wavemax(ft_ratioWave_log2),abs(wavemin(ft_ratioWave_log2)))
	Variable minPVar = wavemin(ft_allTWave)
	minPVar = 10 ^ (floor((log(minPVar))))
	SetAxis/W=$plotName bottom -maxVar,maxVar
	ModifyGraph/W=$plotName mode=3,marker=19,mrkThick=0
	SetDrawEnv/W=$plotName xcoord= bottom,ycoord= left,dash= 3;DelayUpdate
	DrawLine/W=$plotName -maxVar,0.05,maxVar,0.05
	SetDrawEnv/W=$plotName xcoord= bottom,ycoord= left,dash= 3
	DrawLine/W=$plotName 0,1,0,minPVar
	ModifyGraph/W=$plotName zColor(ft_allTWave)={ft_colorWave,0,3,cindexRGB,0,colorTableWave}
	Label/W=$plotName left "P-Value"
	String labelStr = volcanoLabelWave[0] + " / " + volcanoLabelWave[1] + " (Log\\B2\\M)"
	Label/W=$plotName bottom labelStr
	SetWindow $plotName, hook(modified)=filt_thunk_hook
End

Function AddGOTermsToVolcano()
	String plotName = "ftVolcanoPlot"
	WAVE/Z/T filter_GO_Terms
	
	String labelStr
	wfprintf labelStr, "%s\r", filter_GO_Terms
	labelStr = "\Z08" + labelStr
	TextBox/W=$plotName/C/N=GOTs/B=1/F=0/A=LT/X=15.00/Y=0.00 labelStr
End

////////////////////////////////////////////////////////////////////////
// Panel functions
////////////////////////////////////////////////////////////////////////

Function VolcanoIO_Panel()
	WAVE/Z/T volcanoPrefixWave, volcanoLabelWave
	WAVE/Z volcanoParamWave
	if(!WaveExists(volcanoPrefixWave))
		Make/O/N=(2)/T volcanoPrefixWave = {"prefix1*","prefix2*"}
		Make/O/N=(2)/T volcanoLabelWave = {"Test","Control"}
		Make/O/N=(6) volcanoParamWave = {0,0,1,2,1,2}
	endif
	// they are called prefix but suffix (or any wildcard search string is fine)
	String prefix1 = volcanoPrefixWave[0]
	String prefix2 = volcanoPrefixWave[1]
	String label1 = volcanoLabelWave[0]
	String label2 = volcanoLabelWave[1]
	// pick up parameters
	Variable baseVal = volcanoParamWave[0]
	Variable pairOpt = volcanoParamWave[1]
	Variable seedOpt = volcanoParamWave[2]
	Variable foldChange = volcanoParamWave[3]
	Variable seed1 = volcanoParamWave[4]
	Variable seed2 = volcanoParamWave[5]
	
	DoWindow/K VolcanoSetup
	NewPanel/N=VolcanoSetup/K=1/W=(81,73,774,240)
	TitleBox tb0,pos={20,13},size={115,20},title="Search string:",fstyle=1,fsize=11,labelBack=(55000,55000,65000),frame=0
	SetVariable box1,pos={20,33},size={240,16},title="Search string for condition 1 (test):",value=_STR:prefix1
	SetVariable box2,pos={20,54},size={240,16},title="Search string for condition 2 (ctrl):",value=_STR:prefix2
	TitleBox tb1,pos={298,13},size={115,20},title="Labels:",fstyle=1,fsize=11,labelBack=(55000,55000,65000),frame=0
	SetVariable box3,pos={298,33},size={200,16},title="Label for condition 1:",value=_STR:label1
	SetVariable box4,pos={298,54},size={200,16},title="Label for condition 2:",value=_STR:label2
	SetVariable box5,pos={20,115},size={250,16},title="What value represents absent proteins?",format="%g",value=_NUM:baseVal
	SetVariable box6,pos={20,134},size={250,16},title="Fold-change (2 is twofold)?",format="%g",value=_NUM:foldChange
	CheckBox box7,pos={298,131},size={20,20},title="Analyse paired data?",value=pairOpt,mode=0
	CheckBox box8,pos={528,92},size={20,20},title="Reproducibly random?",value=seedOpt,mode=0
	TitleBox tb2,pos={528,13},size={115,20},title="Imputation:",fstyle=1,fsize=11,labelBack=(55000,55000,65000),frame=0
	SetVariable box9,pos={528,33},size={150,16},title="Seed for condition 1:",format="%g",value=_NUM:seed1
	SetVariable box10,pos={528,54},size={150,16},title="Seed for condition 2:",format="%g",value=_NUM:seed2
	TitleBox tb3,pos={528,72},size={115,20},title="Pick a value between 1 and 1000",fstyle=0,fsize=9,frame=0
	Button DoIt,pos={564,140},size={100,20},proc=ButtonProc,title="Do It"
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
	Variable foldChange = V_Value
	ControlInfo/W=$ba.win box7
	Variable pairOpt = V_Value
	ControlInfo/W=$ba.win box8
	Variable seedOpt = V_Value
	ControlInfo/W=$ba.win box9
	Variable seed1 = V_Value
	ControlInfo/W=$ba.win box10
	Variable seed2 = V_Value
	
	Print "Test group:", prefix1, "Labelled", label1, "\rControl group:", prefix2, "Labelled", label2
	Print "Value for imputation:", baseVal, "\rFold-change:", foldChange
	if(pairOpt == 1)
		Print "Using paired data."
	endif
	if(seedOpt == 1)
		Print "Using RNG seeds:", label1, "=", seed1, "&", label2, "=", seed2
	endif
	WAVE/Z/T 	volcanoPrefixWave, volcanoLabelWave
	volcanoPrefixWave[0] = prefix1
	volcanoPrefixWave[1] = prefix2
	volcanoLabelWave[0] = label1
	volcanoLabelWave[1] = label2
	WAVE/Z volcanoParamWave
	volcanoParamWave[0] = baseVal
	volcanoParamWave[1] = pairOpt
	volcanoParamWave[2] = seedOpt
	volcanoParamWave[3] = foldChange
	volcanoParamWave[4] = seed1
	volcanoParamWave[5] = seed2
	VolcanoWorkflowWrapper()
End

Function GOTerm_Panel()
	WAVE/Z/T all_GO_Terms
	WAVE/Z/T all_GO_Terms_sel
	if(!WaveExists(all_GO_Terms) || !WaveExists(all_GO_Terms_Sel))
		return -1
	endif
	String panelName = "GOTPanel"
	KillWindow/Z $panelName
	NewPanel/K=1/N=$panelName/W=(181,179,581,640) as "Select GO Terms"
	ListBox lb1,pos={40,9},size={320,400},listWave=all_GO_terms,selWave=all_GO_Terms_sel,mode=9
	Button reset,pos={40,434},size={100,20},proc=GOTProc,title="Reset"
	Button filter,pos={271,434},size={100,20},proc=GOTProc,title="Filter"
End

Function GOTProc(ba) 
	STRUCT WMButtonAction &ba
	
	if( ba.eventCode != 2 )
		return 0
	endif
	
	WAVE/Z all_GO_Terms_sel
	WAVE/Z/T all_GO_Terms
	
	switch(ba.eventCode)
		case 2:
			if(cmpstr(ba.ctrlName,"reset") == 0)
				all_GO_Terms_sel = 0
				return 0
			elseif(cmpstr(ba.ctrlName,"filter") == 0)
				if(sum(all_GO_Terms_sel) == 0)
					DoAlert 0, "No terms selected"
					return -1
				endif
				// add how many we are filtering for
				Print "Filtering for terms:", sum(all_GO_Terms_sel), "out of", numpnts(all_GO_terms)
				FilteredVolcanoWorkflowWrapper()
				// trigger next part
				return 0
			else
				return -1
			endif
	endswitch
	
	return 0
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
	WAVE/Z volcanoParamWave
	volcanoParamWave[0] = baseVal
	// set param for imputation to random
	volcanoParamWave[2] = 0
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

Function filt_thunk_hook(s)
	Struct WMWinHookStruct& s
	WAVE/Z ft_allTWave
	WAVE/T/Z ft_SHORTNAME
 
	strswitch (s.eventname)
		case "mouseup":

			String s_traceinfo = TraceFromPixel(s.mouseloc.h, s.mouseloc.v, "WINDOW:"+s.winname+";")
			Variable v_pt = str2num(StringByKey("HITPOINT", s_traceinfo))
			String targetTrace = StringByKey("TRACE", s_traceinfo) // now takes whatever trace is clicked on
 
			if (numtype(v_pt) != 2)
				Tag/c/n=t1/b=3/f=0/s=3/v=1/X=10/Y=10 $targetTrace, v_pt, ft_SHORTNAME[v_pt]
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

STATIC Function DeleteEmptyCellsFromTextWave(tw)
	Wave/T tw
	
	Variable numPoints = numpnts(tw)
	
	Variable i
	
	for(i = numPoints - 1; i >= 0; i -= 1)
		if(strlen(tw[i]) == 0)
			DeletePoints i, 1, tw
		endif   
	endfor
End