# TPE-disk
Script for computing stress due to a disk shaped TPE inclusion
%%%%Massimo Nespoli 2022%%%%%%%%%%%%%%%%%%%

The function "TPE_STRESS.m" can be used to compute stresses along a x profile due to a TPE disk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


The "example" directory contains MATLAB scripts to compute the stress field of a TPE inclusion composed of two disks, in a vertical section x,z.
Such an example is the CASE2 reported in:
"Stress changes caused by exsolution of magmatic fluids within an axi-symmetric inclusion" by Maria Elina Belardinelli, Massimo Nespoli and Maurizio Bonafede, GJI (2022)"


Call_double_disk.m : generates the Stress field
TPE_STRESS.m : function called by Call_double_disk.m
Plot_stress.m: Script to plot the computed stress field.
