module PlanetarySystem
using DifferentialEquations
using RecursiveArrayTools

# Newton's constant, sum symbol and number of objects (central star included)
const G = 2.95912208286e-4
const ∑ = sum
const N = 6


# Gravitational newtonian potential
potential(p, t, x, y, z, M) = -G*∑(i->∑(j->(M[i]*M[j])/sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2), 1:i-1), 2:N)


# Solves the N-body problem with a sympletic method, Yoshida6.
# Arguments are the masses, initial velocities and positions and time interval
export NBsolution
function NBsolution(M, vel, pos, tspan)
    nprob = NBodyProblem(potential, M, vel, pos, tspan)
    sol = solve(nprob,Yoshida6(), dt=100)
end


# Write to file


# Read from file


# Nice plot


end # module
