# PlanetarySystem

This Julia package was produced for an assignment of the Computational Physics course held by Erik Schnetter at Perimeter Institute, during the Fall term of 2018.

The package contains tools to solve an N-body problem in Newtonian gravity, and plot and animate the orbits and the consverved quantities. Distance is expressed in A.U., time in Earth days, mass in central star masses.

This package is based on the [DifferentialEquations package tutorial](https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/PhysicalModels/Outer-Solar-System.ipynb).

## Dependencies

Before activating the package, add the following packages:
```
 Pkg.add("DifferentialEquations") 
 Pkg.add("RecursiveArrayTools") 
 Pkg.add("Plots") 
```

## Functions

```
 NBsolution(M::Array{Float64,1}, vel, pos, tspan::Tuple{Float64,Float64})
 ```
To solve the N-Body problem.  vel is the ArrayPartition of the initial velocities of all bodies (see [example](#example)), pos is the ArrayPartition of the initial positions of all bodies, tspan is the time interval.
 
 ```
 myplot(sol, filename::String, planets::Array{String,1})
 ```
To plot the orbit in 3D. sol is the output of NBsolution(), filename is the name of the file where the plot is saved, planets is the array of the planets' names.
 
 ```
 animation(sol, filename::String, planets::Array{String,1})
 ```
 To animate the system in 3D. Expects a .gif file name. Arguments are the same as for the myplot function.
 
 ```
 plot_first_integrals(sol, M::Array{Float64,1}, filename::String, planets::Array{String,1})
 ```
To plot the fractional variation of the conserved quantities (energy and angular momentum vector). All arguments have been described above.

## Example

In this example we solve a planetary system, namely the outer Solar system planets, save a plot of the orbits in position space, save an animation of the system and save a plot of the fractional variation of the conserved quantities. Note that M1>1 takes the inner planets into account.

 ```
using PlanetarySystem
using RecursiveArrayTools

plotfilename="mysol.png"
animfilename="mysol.gif"
fifilename="firstintegrals.png"

M = [1.00000597682, 0.000954786104043, 0.000285583733151, 0.0000437273164546, 0.0000517759138449]
planets = ["Sun", "Jupiter", "Saturn", "Uranus", "Neptune"] 

pos_x = [0.0,-3.5023653,9.0755314,8.3101420,11.4707666]
pos_y = [0.0,-3.8169847,-3.0458353,-16.2901086,-25.7294829]
pos_z = [0.0,-1.5507963,-1.6483708,-7.2521278,-10.8169456]
pos = ArrayPartition(pos_x,pos_y,pos_z)

vel_x = [0.0,0.00565429,0.00168318,0.00354178,0.00288930]
vel_y = [0.0,-0.00412490,0.00483525,0.00137102,0.00114527]
vel_z = [0.0,-0.00190589,0.00192462,0.00055029,0.00039677]
vel = ArrayPartition(vel_x,vel_y,vel_z)

tspan = (0.,200_000.)

sol=NBsolution(M, vel, pos, tspan)
myplot(sol,plotfilename,planets)
animation(sol,animfilename,planets)
plot_first_integrals(sol, M, fifilename, planets)
 ```
 The resulting plots:
 
 ![orbits](https://i.imgur.com/pHHsZg3.png)
 ![orbitsanim](https://i.imgur.com/AfvHrf6.gif)
 ![conservation](https://i.imgur.com/yLKs14W.png)

## Test

The package is tested on random Solar-like systems with a central start initialized at x=(0,0,0) and v=(0,0,0). The test verifies the conservation of energy and angular momentum.

## Acknowledgements 

The author thanks Job Feldbrugge and Stephen Green for their comments and suggestions.

[![Build Status](https://travis-ci.org/laurasberna/PlanetarySystem.jl.svg?branch=master)](https://travis-ci.org/laurasberna/PlanetarySystem.jl)
