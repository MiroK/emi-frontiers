# EMI-FRONTIERS

This repository contains an implementation of the mortar finite element
formulation of the EMI [model](https://doi.org/10.3389/fphy.2017.00048)
needed to reproduce the results of the manufactured test case from Section 3.2
of the paper.

## Dependencies
The following software stack is needed to run the code

- _FEniCS 2017.2.0_ once FEniCS docker infrastructure is in [place](https://docs.docker.com/install/linux/docker-ce/ubuntu/#prerequisites) run
`fenicsproject run quay.io/fenicsproject/stable:2017.2.0`
- _cbc.block_ install from [source](https://github.com/MiroK/cbc.block) in the shell running FEniCS docker image
- _FEniCS_ii_ install from [source](https://github.com/MiroK/fenics_ii) in the shell running FEniCS docker image

## Running the code
There are 2 scripts in the repository. Once in `src` directory

- `python bvp_solver.py` runs a convergence study of the stationary problem
- `python ivp_solver.py` runs a convergence study of the transient problem
