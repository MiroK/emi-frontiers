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
- _tqdm_ install from [source](https://github.com/tqdm/tqdm)

## Running the code
There are 2 scripts in the repository's `src` folder. These run convergence tests
on the intial value problem from Section 3.2 (`mortar_ivp.py`) and a related boundary
value problem (`mortar_bvp.py`). Running `python mortar_bvp.py` the following
errors and rates should be obtained

```
|h               |dim(Wh)         ||ui-uih|_1      ||ue-ueh|_0      ||p-ph|_0        |
|----------------|----------------|----------------|----------------|----------------|
|1.25E-01        | 113            |8.37E-01(nan)   |1.45E+00(nan)   |8.25E-01(nan)   |
|6.25E-02        | 353            |4.32E-01(0.96)  |7.48E-01(0.96)  |2.85E-01(1.53)  |
|3.12E-02        | 1217           |2.18E-01(0.99)  |3.77E-01(0.99)  |9.86E-02(1.53)  |
|1.56E-02        | 4481           |1.09E-01(1.00)  |1.89E-01(1.00)  |3.44E-02(1.52)  |
|7.81E-03        | 17153          |5.45E-02(1.00)  |9.44E-02(1.00)  |1.21E-02(1.51)  |
|3.91E-03        | 67073          |2.73E-02(1.00)  |4.72E-02(1.00)  |4.25E-03(1.51)  |
```
Convergence of the solver for the initial value problem (`python mortar_ivp.py`) is
expected as follows

```
|h               |dim(Wh)         ||u-uh|_1        ||u-uh|_0        ||p-ph|_0        ||u-uh|_oo       ||duh-du|_oo     |
|----------------|----------------|----------------|----------------|----------------|----------------|----------------|
|1.25E-01        | 113            |2.12E+00(nan)   |1.03E-01(nan)   |4.56E-01(nan)   |1.80E-01(nan)   |1.23E-01(nan)   |
|6.25E-02        | 353            |1.11E+00(0.94)  |2.81E-02(1.88)  |1.45E-01(1.66)  |6.14E-02(1.55)  |4.81E-02(1.35)  |
|3.12E-02        | 1217           |5.59E-01(0.98)  |7.20E-03(1.97)  |4.00E-02(1.85)  |1.98E-02(1.63)  |1.70E-02(1.50)  |
|1.56E-02        | 4481           |2.80E-01(1.00)  |1.81E-03(1.99)  |1.05E-02(1.94)  |6.16E-03(1.68)  |5.64E-03(1.59)  |
|7.81E-03        | 17153          |1.40E-01(1.00)  |4.54E-04(2.00)  |2.66E-03(1.98)  |1.86E-03(1.73)  |1.79E-03(1.66)  |
|3.91E-03        | 67073          |7.02E-02(1.00)  |1.14E-04(2.00)  |6.68E-04(1.99)  |5.48E-04(1.76)  |5.44E-04(1.71)  |
|1.95E-03        | 265217         |3.51E-02(1.00)  |2.84E-05(2.00)  |1.67E-04(2.00)  |1.58E-04(1.79)  |1.61E-04(1.76)  |
```
