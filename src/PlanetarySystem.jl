module PlanetarySystem
using DifferentialEquations
using RecursiveArrayTools
#using Plots
#using HDF5
#using ImageMagick

# Newton's constant, sum symbol and number of objects (central star included)
const G = 2.95912208286e-4
const ∑ = sum
const N = 5#6


# Gravitational newtonian potential
potential(p, t, x, y, z, M) = -G*∑(i->∑(j->(M[i]*M[j])/sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2), 1:i-1), 2:N)


# Solves the N-body problem with a sympletic method, Yoshida6.
# Arguments are the masses, initial velocities and positions and time interval
export NBsolution
function NBsolution(M, vel, pos, tspan)
    nprob = NBodyProblem(potential, M, vel, pos, tspan)
    sol = solve(nprob,Yoshida6(), dt=100)
end


#=Write to file in HDF5 form
export writefile
function writefile(sol, filename::String)
    h5write(filename, "data", sol)
    end

#Read back
export readfile
function readfile(filename::String)
    data=h5read(filename, "data")
end=#


#=Plot
export myplot
function myplot(sol, planets::Array{String,1}, filename::String, mytitle::String)
    #pyplot()
    orbitplot(sol,
              body_names=planets,
              xlim=(-0.01,0.01), ylim=(-0.01,0.01), zlim=(-0.005,0.005),
              w=1.5,
              xlabel = "x", ylabel = "y", zlabel = "z",
              legend=:bottomleft,
              size=(800, 600),
              title  = mytitle)
    savefig(filename)
    end =#

#Animation
export animation
function animation(sol, filename::String)
    animate(sol,lw=3,every=100,filename)
end

end # module
