# The package contains tools to solve an N-body problem in Newtonian gravity,
# and plot and animate the orbits and the consverved quantities.
# This package is based on the DifferentialEquations package tutorial:
# https://github.com/JuliaDiffEq/DiffEqTutorials.jl/blob/master/PhysicalModels/Outer-Solar-System.ipynb

module PlanetarySystem

using DifferentialEquations
using RecursiveArrayTools
using Plots


const G = 2.95912208286e-4 # Newton's constant, units: [L]= AU, [M]= M_sun, [t]= Earth days 
const ∑ = sum # symbol for sum function
const N = 5 #Number of bodies (central star included)


#############################################################

# Physical model:

# Gravitational newtonian potential for the NBsolution() solver.
# p is a parameter, t the time, the coordinates x, y, z are arrays running over the bodies,
# M is the array of the bodies' masses
potential(p, t, x, y, z, M) = -G*∑(i->∑(j->(M[i]*M[j])/sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2), 1:i-1), 2:N)

# Kinetic energy
# the velocity components  vx, vy, vz are 2D arrays running over the bodies and the time steps,
# M is the array of the bodies' masses
kinetic(vx, vy, vz, M) = 1/2 * ( ∑(M.*vx.*vx,dims=1) + ∑(M.*vy.*vy,dims=1) + ∑(M.*vz.*vz,dims=1) )
# potential, to compute the Hamiltonian
# the coordinates x, y, z are 2D arrays running over the bodies and the time steps,
# M is the array of the bodies' masses
pot(x, y, z, M)= - G*∑(i->∑(j->(M[i]*M[j])./sqrt.((x[i,:]-x[j,:]).^2+(y[i,:]-y[j,:]).^2+(z[i,:]-z[j,:]).^2), 1:i-1), 2:N)

# Hamiltonian (total energy) - variables as above
H(x, y, z, vx, vy, vz, M) = kinetic(vx, vy, vz, M) + pot(x, y, z, M)'

# Angular momentum components  - variables as above
Lx(x, y, z, vx, vy, vz, M) =  ∑(i-> M[i] .*  (vy[i,:] .* x[i,:] - vx[i,:] .* y[i,:]), 1:N)'
Ly(x, y, z, vx, vy, vz, M) =  ∑(i-> M[i] .*  (-vz[i,:] .* x[i,:] + vx[i,:] .* z[i,:]), 1:N)'
Lz(x, y, z, vx, vy, vz, M) =  ∑(i-> M[i] .*  (vz[i,:] .* y[i,:] - vy[i,:] .* z[i,:]), 1:N)'  

export H, Lx, Ly, Lz

#############################################################

# Main
export SolveAndPlot
function SolveAndPlot(M::Array{Float64,1}, vel, pos, tspan::Tuple{Float64,Float64}, planets::Array{String,1}, plotfilename::String, animfilename::String, fifilename::String)

sol=NBsolution(M, vel, pos, tspan)
myplot(sol,plotfilename,planets)
animation(sol,animfilename,planets)
plot_first_integrals(sol, M, fifilename, planets)

end


#############################################################


# Solves the N-body problem with a sympletic method, Yoshida6.

# M is the array of the bodies' masses,
# vel contains the initial velocity vectors of all bodies (see example in Readme)
# pos contains the initial position of all bodies (see example in Readme)
# tspan is the time interval

export NBsolution
function NBsolution(M::Array{Float64,1}, vel, pos, tspan::Tuple{Float64,Float64})
    nprob = NBodyProblem(potential, M, vel, pos, tspan)
    sol = solve(nprob, Yoshida6(), dt=10) # can change the method or the time step dt
end

#############################################################

# Plots the orbits in coordinate space and saves to file
# sol is the output of NBsolution()
# filename is the name of the file where the plot is saved
# planets is the array of the planets' names

export myplot
function myplot(sol, filename::String, planets::Array{String,1})
    plot(size=(500, 400),
         title  = "Planetary orbits",
         legend=:bottomleft,
         xlim=(-30,30), ylim=(-30,30), zlim=(-30,30)
         ,markerstrokewidth=0
         #,bg=:black
         #,xlabel = "x", ylabel = "y", zlabel = "z"
         )
    
    for i in 1:N
    scatter!(sol,
              vars=(3*N+i,4*N+i,5*N+i),
              label=planets[i],
              markersize=1,
              markerstrokewidth=0)
    end
          
    savefig(filename)
end

#############################################################

# Animates the planets and saves to file
# see above

export animation
function animation(sol, filename::String, planets::Array{String,1})
    
    nplots=500 
    element(i)=round(Int, length(sol.t)/nplots-1)*i
    
    anim = @animate for i=1:nplots #    
        scatter((sol[3*N+1,element(i)],sol[4*N+1,element(i)],sol[5*N+1,element(i)]),
                #xlabel = "x", ylabel = "y", zlabel = "z",
                xlim=(-30,30), ylim=(-30,30), zlim=(-30,30),
                label=planets[1],
                markerstrokewidth=0,
                markercolor=:yellow,
                marker=10,
                title  = "Planetary system",
                size=(500, 400))
        if(N>1)
            for j in 2:(N)
        scatter!((sol[3*N+j,element(i)],sol[4*N+j,element(i)],sol[5*N+j,element(i)]),
                 marker = 5,
                 markerstrokewidth=0,
                 label=planets[j])
            end
        end
    end
   
    gif(anim, filename, fps = 17) #change fps to adjust the speed of the animation

end

#############################################################

# Plots the fractional variation of conserved quantities (energy and angular momentum)
# see above

export plot_first_integrals
function plot_first_integrals(sol, M::Array{Float64,1}, filename::String, planets::Array{String,1})
    
    vx=sol[1:N,:]
    vy=sol[N+1:2*N,:]
    vz=sol[2*N+1:3*N,:]
    x= sol[3*N+1:4*N,:]
    y= sol[4*N+1:5*N,:]
    z= sol[5*N+1:6*N,:]
    
    plot(sol.t, (H(x, y, z, vx, vy, vz, M)[1] .- H(x, y, z, vx, vy, vz, M)')./H(x, y, z, vx, vy, vz, M)[1],
         lab="Energy variation",
         title="Conserved Quantities",
         xlabel="t",
         size=(500, 400)) 
    plot!(sol.t, (Lx(x, y, z, vx, vy, vz, M)[1].-Lx(x, y, z, vx, vy, vz, M)')./Lx(x, y, z, vx, vy, vz, M)[1], lab="Angular momentum variation, x")
    plot!(sol.t, (Ly(x, y, z, vx, vy, vz, M)[1].-Ly(x, y, z, vx, vy, vz, M)')./Ly(x, y, z, vx, vy, vz, M)[1], lab="Angular momentum variation, y")
    plot!(sol.t, (Lz(x, y, z, vx, vy, vz, M)[1].-Lz(x, y, z, vx, vy, vz, M)')./Lz(x, y, z, vx, vy, vz, M)[1], lab="Angular momentum variation, z")

    savefig(filename)
   
end



end # module
