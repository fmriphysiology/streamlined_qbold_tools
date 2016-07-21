streamlined_qbold_tools:

Scripts used to produce quantitative maps of baseline brain oxygenation from FLAIR-GASE data. See streamlined-qBOLD paper for more info.

NOTES:

run.qbold_make is a wrapper to run the example data ... email me for this if your interested 
Add x.zaverage and x.qbold_make to shell path
Add path to Matlab (x.qbold_make, line 48) so you can run from command line
Add the ase_qbold_3d function to matlab path in setup.m
Possible changes ...
    - give option for output directory name in x.qbold_make
    - don't bother with the number counter when looping over voxels in ase_qbold_3d
    - upload example image?
