# PlanetarySystem

This Julia package was produced for an assignment of the Computational Physics course held by Erik Schnetter at Perimeter Institute, during the Fall term of 2018.

The package solves the N-body problem in Newtonian gravity and is based on the [DifferentialEquations package tutorial](https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/PhysicalModels/Outer-Solar-System.ipynb).

## Dependencies

```
 Pkg.add("DifferentialEquations") 
 Pkg.add("RecursiveArrayTools") 
 Pkg.add("Plots") 
 Pkg.add("LinearAlgebra")
```

## Functions

```
 NBsolution(M::Array{Float64,1}, vel, pos, tspan::Tuple{Float64,Float64})
 ```
 To solve the N-Body problem.
 
 ```
 myplot(sol, filename::String, planets::Array{String,1}, mytitle::String)
 ```
 To plot the orbit in 3D.
 
 ```
 animation(sol, filename::String, planets::Array{String,1}, mytitle::String)`
 ```
 To animate the system in 3D. Expects a .gif file name.
 
 ```
 plot_first_integrals(sol, M::Array{Float64,1}, filename::String, planets::Array{String,1})
 ```
To plot the fractional variation of the conserved quantities (energy and angular momentum vector).

## Acknowledgements 

The author thanks Job Feldbrugge and Stephen Green for their comments and suggestions.

[![Build Status](https://travis-ci.org/laurasberna/PlanetarySystem.jl.svg?branch=master)](https://travis-ci.org/laurasberna/PlanetarySystem.jl)
