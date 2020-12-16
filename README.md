# Modematching
Phonon normal mode-matching algorithm used for upscaling low-accuracy semi-empirical (DFTB) phonons to full ab initio (DFT)


##
This code contains a full suite of finite temperature free energy predictions based on pre-calculated harmonic phonons. 

##
The base functionality takes the output of a Gamma-point DFT and DFTB calculation (ref.yaml, shift.yaml files). These are required to contain the normalized eigenvectors and frequencies, and will serve as the initial matching assignnment. A third (ss.yaml) file containing the frequencies of the raw DFTB supercell calculation containing the entire density of states is required. The differences in the frequencies of the 2x G-pt calculations will be added here as the "shift", and the full integration of the shifted density of states (DOS) will serve as the upscaled free energy.

##
An additional correction to the DOS can be applied with the inclusion of acoustic modes via elastic constants and toggled by the $/ACOUSTICS card. This will require either an input of pre-calculated Elastic Constants (.ECs file) or stress/strain files for the ECs to be evaluated here (Toggle CALC ECs = TRUE in $/Acoustics, files .stress and .strain). This will ensure convergence of the acoustic modes near Gamma. 

##
This program will also allow for the inclusion of multiple sets of frequencies for Quasiharmonic Approximations (QHA), predicting accurate thermal expansions and sublimation enthalpies. Simply put as many frequency files in the working directory as you wish to be included in the Fvib(V) calculation, as well as the volumes in the same order in the input file, and directories will be made for each match/acoustic treatment. 

##
Finally, ModeMatch can predict solid-state phase diagrams given multiple structures in the working directory. Simply identify multiple polymorphs by listing more structures in the infile inder $/Structures, and name the input files accordingly. 

### Running the Code
The code requires 2 things in the working directory:
1. Input file titled `infile` containing job control
2. Directory titled `datafiles` containing all of the necessary data files for the job you're running. More on this in the **Input File** section. 

### Input File
