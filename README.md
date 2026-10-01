
~~~text
███    ███ ██   ██ ██ ████████ ██████   ██████  
████  ████ ██   ██ ██    ██         ██ ██       
██ ████ ██ ███████ ██    ██     █████  ███████  - Fibra
██  ██  ██ ██   ██ ██    ██         ██ ██    ██ 
██      ██ ██   ██ ██    ██    ██████   ██████         
~~~


#### GPU-based Finite difference code for DNS of Multiphase Homogenous isotropic turbulence with a single fibre

Developers:
* A. Roccon (MHIT36 + porting of fiber tracking from FluTAS to MHIT36)
* V. Agrawal (Original implementation of Fiber tracking in FluTAS)

Tested on:
* Milton (1 x RTX5000)
* Leonardo (1 x A100)
* Marge (1 X RTX6000 QMax)

Current capabiltiies:
* DNS of single-phase flow
* Tracking of a stiff fiber using the method of Agrawal et al. 2024, DOI: https://doi.org/10.1016/j.cma.2023.116495
* Extenension to many fibers (planned, TBD)

#### Systems supported:
* Unix + nvfortran 

#### To run a simulation:
* go to src folder and run ./compile_local.sh (if you have one GPU, UNIX system) or ./compile_Leo.sh to run on Leondaro supercomputer.

#### Parallelization strategy
* The code is serial and exploit a single GPU (GPU-resident)

#### Output files.
* Files containing the Eulerian fields (u\_\*\*\*, v\_\*\*\*\*, w\_\*\*\*\*, p\_\*\*\*\* and phi\_\*\*\*\*) and the fibre positions and orientations (fib\_\*\*\*) are stored in src/output


### Structural solver validation
## Free-free bending frequencies vs Euler–Bernoulli theory

Fiber: L = 1, d = 0.025, ρ_f = 10, E = 218·10³, G = 83.846·10³.
Settings common to all runs (`fib_param.f90`): `fib_coupling = .false.`, `g_vec = 0`,
`pert_amp = 1.d-3`, `fib_log = 1`. `input.inp`: dt = 2.5·10⁻⁴, 8000 steps.

Analytical: ω_n = (β_n L)² √(EI / (ρ_f A L⁴)), β_n L = 4.7300, 7.8532, 10.9956 (modes 1–3).
With `align_frame = .false.` (FluTAS frame) bending in the y–z plane uses G·J instead of E·I.

| Case | Mode | `align_frame` | `pert_dir` | `nel` | Stiffness | ω analytical | ω GPU | Error GPU | ω CPU (reference) | `fiber_freq.py` args |
|---|---|---|---|---|---|---|---|---|---|---|
| 1  | 1 | `.false.` | (0,1,0) | 12 | G·J | 18.108  | 17.903  | −1.13% | 17.903  | `2.5e-4 1 1` |
| 2  | 1 | `.false.` | (1,0,0) | 12 | E·I | 20.646  |         |        | 20.200  | `2.5e-4 0 1` |
| 3  | 2 | `.false.` | (0,1,0) | 12 | G·J | 49.915  |         |        | 50.562  | `2.5e-4 1 2` |
| 4  | 3 | `.false.` | (0,1,0) | 12 | G·J | 97.853  |         |        | 85.070  | `2.5e-4 1 3` |
| 5  | 1 | `.true.`  | (0,1,0) | 12 | E·I | 20.646  |         |        | 20.189  | `2.5e-4 1 1` |
| 6  | 1 | `.true.`  | (0,1,0) | 24 | E·I | 20.646  |         |        | 20.565  | `2.5e-4 1 1` |
| 7  | 1 | `.true.`  | (0,1,0) | 48 | E·I | 20.646  |         |        | 20.636  | `2.5e-4 1 1` |
| 8  | 2 | `.true.`  | (0,1,0) | 12 | E·I | 56.912  |         |        | 57.677  | `2.5e-4 1 2` |
| 9  | 2 | `.true.`  | (0,1,0) | 24 | E·I | 56.912  |         |        | 56.844  | `2.5e-4 1 2` |
| 10 | 2 | `.true.`  | (0,1,0) | 48 | E·I | 56.912  |         |        | 56.700  | `2.5e-4 1 2` |
| 11 | 3 | `.true.`  | (0,1,0) | 12 | E·I | 111.570 |         |        | 83.971  | `2.5e-4 1 3` |
| 12 | 3 | `.true.`  | (0,1,0) | 24 | E·I | 111.570 | 110.245 | −1.19% | 110.245 | `2.5e-4 1 3` |
| 13 | 3 | `.true.`  | (0,1,0) | 48 | E·I | 111.570 |         |        | 110.430 | `2.5e-4 1 3` |

