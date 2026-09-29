
~~~text
███    ███ ██   ██ ██ ████████ ██████   ██████  
████  ████ ██   ██ ██    ██         ██ ██       
██ ████ ██ ███████ ██    ██     █████  ███████  - Fibre
██  ██  ██ ██   ██ ██    ██         ██ ██    ██ 
██      ██ ██   ██ ██    ██    ██████   ██████         
~~~


#### GPU-based Finite difference code for DNS of Multiphase Homogenous isotropic turbulence (Single fibre)

Main developer: A. Roccon 

Tested on:
* Milton (1 x RTX5000)
* Leonardo (1 x A100)
* Marge (1 X RTX6000 QMax)

Current capabiltiies:
* DNS of single-phase flow
* Tracking of a stiff fiber using the method of Agrawal et al. 2024, DOI: https://doi.org/10.1016/j.cma.2023.116495

#### Systems supported:
* Unix + nvfortran 

#### To run a simulation:
* go to src folder and run ./compile_local.sh (if you have one GPU, UNIX system) or ./compile_Leo.sh to run on Leondaro supercomputer.

#### Parallelization strategy
* The code is serial and exploit a single GPU (GPU-resident)

#### Output files.
* Files containing the Eulerian fields (u\_\*\*\*, v\_\*\*\*\*, w\_\*\*\*\*, p\_\*\*\*\* and phi\_\*\*\*\*) and the particle positions (xp\_\*\*\*) are stored in src/output

