# PlanetarySystem

This Julia package was produced for an assignment of the Computational Physics course held by Erik Schnetter at Perimeter Institute, during the Fall term of 2018.

The package solves the N-body problem in Newtonian gravity and is based on the [DifferentialEquations package tutorial](https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/PhysicalModels/Outer-Solar-System.ipynb).

## Dependencies

```
 Pkg.add("DifferentialEquations"), Pkg.add("RecursiveArrayTools"), Pkg.add("Plots"), Pkg.add("LinearAlgebra")
```

## Functions

```
 NBsolution(M, vel, pos, tspan)
 ```

## Acknowledgements 

The author thanks Job Feldbrugge and Stephen Green for their comments and suggestions.

[![Build Status](https://travis-ci.org/laurasberna/PlanetarySystem.jl.svg?branch=master)](https://travis-ci.org/laurasberna/PlanetarySystem.jl)
