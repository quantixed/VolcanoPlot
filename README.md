# VolcanoPlot
Making volcano plots from proteomic data in IgorPro

Pairwise comparison of "intensity" data from proteomics data from two conditions (MaxQuant).

Run from Macros menu and tell Igor which two conditions you want to compare (give Igor a prefix for waves), e.g. `ctrl*` and `test*` will find ctrl_1,ctrl_2 etc.

Also tell Igor what the basevalue is. This is the intensity assigned to any proteins which were not detected. NaN is currently not supported.

The protein names need to be in TextWaves with names `NAME` and `SHORTNAME`. I will add more flexibility to this soon.

Igor will make a Volcano Plot of this comparison. Proteins are coloured according to magnitude of change from condition 2 (2^2) and according to P-value. Proteins can be clicked on to reveal their `SHORTNAME`.

Minimum data to run (all in `root`):

- Three waves of condition1
- Three waves of condition2 - logical naming advised
- NAME and SHORTNAME textwaves