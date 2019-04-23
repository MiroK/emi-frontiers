# EMI-FRONTIERS

This repository contains an implementation of the mortar finite element
formulation of the EMI [model](https://doi.org/10.3389/fphy.2017.00048)
needed to reproduce the results of the manufactured test case from Section 3.2
of the paper.

## Dependencies
The following software stack is needed to run the code

- FEniCS 2017.2.0: once FEniCS docker infrastructure is in [place](fenicsproject run quay.io/fenicsproject/stable:2017.2.0) run
`fenicsproject run quay.io/fenicsproject/stable:2017.2.0`
- cbc.block: install from [source](https://github.com/MiroK/cbc.block) in the shell running FEniCS docker image
- [FEniCS_ii]: install from [source](https://github.com/MiroK/fenics_ii) in the shell running FEniCS docker image

## Running the code
There are 2 scripts in the repository
- `python bwe_problem.py` runs a convergence study of the stationary problem
- `python ivp_problem.py` runs a convergence study of the transient problem