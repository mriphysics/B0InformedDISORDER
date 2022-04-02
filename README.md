# B0InformedDISORDER

This repository provides the tools to implement the methods and reproduce the experiments presented in the manuscript 'Data-driven motion-corrected brain MRI incorporating pose dependent B0 fields', Y Brackenier, L Cordero-Grande, T Tomi-Tricot, T Wilkinson, P Bridgen, AN Price, SJ Malik, E De Vita and JV Hajnal, 2022.

The code has been developed in MATLAB and has the following structure:

###### ./Control
Contains the scripts that store parameters for setting the paths and generating the figures.

###### ./Data
contains the data needed to reproduce all experiments.

###### ./Data/Exp1_Pose
contains the Nifti files of the registered pose-acquisitions at 3T and 7T. A prepared reconstruction structure is included which is used in the simulations.

###### ./Data/Exp2_Motion
contains the prepared reconstruction structures of the controlled motion experiments.

###### ./Manuscript
contains the scripts to generate all the figures published in the manuscript.

###### ./Methods
contains all the functions needed to perform motion and/or B0 correction.

###### ./Methods/DISORDER
contains the original DISORDER repository as presented in https://github.com/mriphysics/DISORDER. 

###### ./Methods/DISORDER_Modified
contains the functions from https://github.com/mriphysics/DISORDER that have been modified in the current reconstruction framework. Corresponding functions have been deleted in ./Methods/DISORDER 

###### ./Methods/DISORDER_Added
contains additional implemented functions needed for the proposed reconstruction.

###### ./Methods/Utilities
contains a set of helper functions used throughout the repository.

###### ./Results
contains the results from running the scripts in ./Scripts. Results for the pose experiment, motion experiment and simulations are stored in respectively ./Results/Exp1_Pose, ./Results/Exp2_Motion and ./Results/Simulations.
For the in-vivo experiments, Nifti files are generated in subfolder *An-Ve* with suffixes **_Aq_MotCorr.nii* (no motion correction), **_Di_MotCorr.nii* (motion correction) and **_Di_MotB0Corr.nii* (motion + B0 correction).

###### ./Scripts
contains the script that will run experiment 1, experiment 2 and the simulations.

NOTE 1: For simulations, results are stored in Matlab structures (and not Nifti images). These structures will be read in when generating the figures.

NOTE 2: Before running the scripts in ./Manuscript, all scripts in ./Scripts should be executed to generate the necessary results.

NOTE 3: Due to the current ethics regulations, the data required to reproduce the results is not share online. Please contact yannick.brackenier@kcl.ac.uk if you are interested in the data. 

