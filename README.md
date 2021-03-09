# Modematching
Phonon normal mode-matching algorithm used for upscaling low-accuracy semi-empirical (DFTB) phonons to full ab initio (DFT)

<p align="center">
  <img src="https://github.com/cjcook41/Modematching/blob/cjcook41-patch-3/images/GitFig1.png">
</p>

This code contains a full suite of finite temperature free energy predictions based on pre-calculated harmonic phonons. 


The base functionality takes the output of a Gamma-point DFT and DFTB calculation (`xtal.ref.yaml`, `xtal.shift.yaml` files). These are required to contain the normalized eigenvectors and frequencies, and will serve as the initial matching assignnment. A third (`xtal.ss.yaml`) file containing the frequencies of the raw DFTB supercell calculation containing the entire density of states is required. The differences in the frequencies of the two Gamma-pt calculations will be added here as the "shift", and the full integration of the shifted density of states (DOS) will serve as the upscaled free energy.


An additional correction to the DOS can be applied with the inclusion of acoustic modes via elastic constants and toggled by the `\$ACOUSTICS` card. This will require either an input of pre-calculated Elastic Constants (`xtal.ecs` file) or stress/strain files for the ECs to be evaluated here (Toggle `CALC ECs  TRUE` in `\$Acoustics`, files `xtal.stress` and `xtal.strain`). This will ensure convergence of the acoustic modes near Gamma. 

This program will also allow for the inclusion of multiple sets of frequencies for Quasiharmonic Approximations (QHA), predicting accurate thermal expansions and sublimation enthalpies. Simply put as many frequency files in the working directory as you wish to be included in the Fvib(V) calculation (For n points: xtal1 - xtaln), as well as the volumes in the same order in the input file, and directories will be made for each match/acoustic treatment. 

Finally, ModeMatch can predict solid-state phase diagrams given multiple structures in the working directory. Identify multiple polymorphs by listing more structures in the infile inder `\$crystals`, and name the input files accordingly.

## Installation

### Dependencies
* numpy (>= 1.18.1)
* sympy (>= 1.5.1)
* pandas (>=1.0.1)

### From pip
First, it is recommended that a virtual environment is created for module control, e.g.:

`conda create -n modematch`

Modematch can then be installed via pip:

`conda activate modematch`

`pip install modematch`


### From source

`git clone https://github.com/cjcook41/Modematching.git`

`conda activate modematch`

`python setup.py sdist`

`pip install modematch/modematch-0.0.4.tar.gz/`


## Running the Code
The code requires 2 things in the working directory:
1. Input file titled `infile` containing job control
2. Directory titled `datafiles` containing all of the necessary data files for the job you're running. More on this in the **Input File** section. 

Run the following commands from a python3 interpreter:
* `from modematch import main`
* `main.RunJob()`

### Elastic Properties
Control of elastic calcualations is listed under the `\$acoustics` card. `Christoffel true` will correct raw frequencies using elastic constants and Christoffel matrices. `EC true` will evaluate the elastic constants on the fly using input stresses and strains.  

### Phase Diagrams
For phase diagrams, add another polymorph name under the `\$crystals` card. **NOTE**: The names listed under here will need to be the prefix of the input data files! `\$P2atoms` and `\$P2volumes` will need to be specified as well, containing atomic information and cell volumes of the second polymorph specified. 
Volume dependent free energies are evaluated at a series of pressures up to `\$Pmax` in GigaPascals (0.01GPa steps),  and temperatures (0 - 500K, 10K steps), and third order polynomial fits are used to determine the transition temperature at some pressure (where relative thermodynamic stability inverts). "No tranition" will result in imaginary temperatures. 

### Sample Input File
`infile` is a sample input file for Phase I carbon dioxide. This input file will perform the base normal mode shift, evaluate free energies at multiple volumes of an Alpha `xtal`, evaluate the corrected acoustic modes using pre-calculated elastic constants, and finally perform the quasiharmonic approximation for thermal expansion. 


### Outputs
* For 1 polymorph, 1 volume (No QHA):
  * Thermodynamic data / absolulte free energies

* 1 polymorph, multiple volumes (QHA):
  * `xtal.csv` file containing relevent free energies and predicted cell volumes as a function of T

* 2 polymorphs, multiple volumes (QHA + Phase Diagram):
  * `xtal.csv` file containing relevent free energies and predicted cell volumes as a function of T
  * PhaseDiagram file containing temperatures, enthalpies, and entropies of transition as a function of P
  
  
## Citation
Please refer to the following publications for more information and citation:
* https://doi.org/10.1063/5.0032649
