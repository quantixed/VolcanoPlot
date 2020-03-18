# VolcanoPlot
Making volcano plots from proteomic data in IgorPro

## From MaxQuant

You can load data from a `proteinGroups.txt` file using the *Load MaxQuant Data...* option. This will load the data and do the analysis.
This procedure will do a pairwise comparison of "LFQ Intensity" from proteomics data from two conditions.

Tell Igor which two conditions you want to compare (give Igor a prefix for waves), e.g. `ctrl*` and `test*` will find ctrl_1,ctrl_2 etc. and test_1,test_2 etc.

Also tell Igor what the basevalue is. This is the intensity assigned to any proteins which were not detected. Default is 0 (for MaxQuant).

## Manual analysis

Analysis can also be started manually from the Macros menu.
In this case, supply data (intensity, LFQ_Intensty or peptides) for the conditions you want to analyse (with logical naming) and also the protein names. These need to be in TextWaves with names `NAME` and `SHORTNAME`.

Minimum data to run (all in `root:`):

- Three waves of condition1
- Three waves of condition2 - logical naming advised
- NAME and SHORTNAME textwaves

## Outputs

Igor will make a Volcano Plot of the comparison. Proteins are coloured according to magnitude of change from control and according to P-value. Proteins can be clicked on to reveal their `SHORTNAME`.

Transforms and imputation are done exactly as described for the default settings in Perseus.

## Further notes

It is possible to do pairwise comparisons for the volcano plot by clicking the checkbox.

A PCA is generated for the selected data. If you want to compare some other data using PCA use the menu item. This is useful if you have more than two conditions and want to do a PCA on everything.