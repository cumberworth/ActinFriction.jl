# Actin rings dynamics simulation package

A Julia package for simulating the dynamics of passively crosslinked actin rings.

This package allows the SDEs described in [Ref. 1](#references) to be solved, and provides methods for directly calculating the friction coefficients described in the same paper.

## Installation

The package can be installed by starting the Julia REPL, typing `]` to enter package mode, and running
```
add ActinFriction
```
to install from the General registry, or by running
```
add https://github.com/cumberworth/ActinFriction.jl
```
to install directly from the development repository.

A related python package, [actinfrictionpy](https://github.com/cumberworth/actinfrictionpy), includes code for analyzing and plotting the output from these ...

## References

[1] A. Cumberworth and P. R. ten Wolde, Constriction of actin rings by passive crosslinkers, [arXiv:2203.04260 [physics.bio-ph]](https://doi.org/10.48550/arXiv.2203.04260).

## Links

[Julia programming language](https://julialang.org/)

[Plotting package](https://github.com/cumberworth/actinfrictionpy)

[Replication package Ref. 1](https://doi.org/10.5281/zenodo.6327217)
