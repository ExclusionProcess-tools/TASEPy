# TASEPy
TASEPy is a python package providing a numerical solution for the inhomogeneous Totally Asymmetric Simple Exclusion Process (TASEP) in the stationary state. The TASEP is a paradigmatic lattice model for one-dimensional particle transport subject to excluded-volume interactions. The package TASEPy provides an implementation of the Power Series Approximation (PSA) introduced in [1,2] for solving the TASEP in the limit in which the initiation rate at which new particles are added to the lattice is small. The iterative solution proposed here is detailed in the affiliated paper [3,4].

All numerical methods are implemented in Python 3.

## Tutorial

A short tutorial demonstrating TASEPy is given in this [python notebook](<tutorial_TASEPy.ipynb>).

## Benchmarks

The TASEPy has been tested using exact results obtained by solving the master equation for small systems, and numerical results obtained by stochastic simulations for large systems. These tests can be found in this [python notebook](<benchmarks_TASEPy.ipynb>).

## Bibliography
[1] J. Szavits-Nossan, L. Ciandrini and M. C. Romano, Deciphering mRNA Sequence Determinants of Protein Production Rate, *Physical Review Letters* 120, 128101 (2018) \[[url](https://doi.org/10.1103/PhysRevLett.120.128101)\]

[2] J. Szavits-Nossan, M. C. Romano and L. Ciandrini, Power series solution of the inhomogeneous exclusion process, *Physical Review E* 97, 052139 (2018) \[[url](https://doi.org/10.1103/PhysRevE.97.052139)\]

[3] L. Ciandrini, R. Crisostomo, J. Szavits-Nossan, TASEPy: a Python-based package to iteratively solve the inhomogeneous exclusion process, *SciPost Physics Codebases, 022*  \[[url](https://www.scipost.org/SciPostPhysCodeb.22)\] (2023)

[4] L. Ciandrini, R. Crisostomo, J. Szavits-Nossan, Codebase release 1.1 for TASEPy, *SciPost Physics Codebases, 022*  \[[url](https://www.scipost.org/SciPostPhysCodeb.22-r1.1)\] (2023)

