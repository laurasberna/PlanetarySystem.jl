module PlanetarySystem
using DifferentialEquations
using RecursiveArrayTools
using Plots
using LinearAlgebra
#using HDF5
#using ImageMagick

# Newton's constant, sum symbol and number of objects (central star included)
const G = 2.95912208286e-4
const ∑ = sum
const N = 5#6


# Gravitational newtonian potential
potential(p, t, x, y, z, M) = -G*∑(i->∑(j->(M[i]*M[j])/sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2), 1:i-1), 2:N)
#Kinetic energy
kinetic(vx, vy, vz, M) = 1/2 * ( ∑(M.*vx.*vx,dims=1) + ∑(M.*vy.*vy,dims=1) + ∑(M.*vz.*vz,dims=1) )
#potential again
pot(x, y, z, M)= - G*∑(i->∑(j->(M[i]*M[j])./sqrt.((x[i,:]-x[j,:]).^2+(y[i,:]-y[j,:]).^2+(z[i,:]-z[j,:]).^2), 1:i-1), 2:N)                                         
#Hamiltonian
export H
H(x, y, z, vx, vy, vz, M) = kinetic(vx, vy, vz, M) + pot(x, y, z, M)'
#Angular momentum
L(x, y, z, vx, vy, vz, M) = ∑(i-> M[i] .* sqrt.( (vy[i,:] .* x[i,:] - vx[i,:] .* y[i,:]).^2 + (-vz[i,:] .* x[i,:] + vx[i,:] .* z[i,:]).^2 + (vz[i,:] .* y[i,:] - vy[i,:] .* z[i,:]).^2 ), 1:N)'

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


#Plot of the orbits in coordinate space
export myplot
function myplot(sol, filename::String, planets::Array{String,1}, mytitle::String)
    plot(size=(800, 600),
         title  = mytitle,
         legend=:bottomleft,
         xlim=(-30,30), ylim=(-30,30), zlim=(-30,30),
         bg=:black
         #,xlabel = "x", ylabel = "y", zlabel = "z"
         )
    
    for i in 1:N
    plot!(sol,
              vars=(3*N+i,4*N+i,5*N+i),
              label=planets[i],
              w=1)
    end
          
    savefig(filename)
    end 

#Animation
export animation
function animation(sol, filename::String, planets::Array{String,1}, mytitle::String)
    gr()
    element(i)=4*i
    
    anim = @animate for i=1:500
        
        scatter((sol[3*N+1,element(i)],sol[4*N+1,element(i)],sol[5*N+1,element(i)]),
                #xlabel = "x", ylabel = "y", zlabel = "z",
                xlim=(-30,30), ylim=(-30,30), zlim=(-30,30),
                label=planets[1],
                markerstrokewidth=0,
                markercolor=:yellow,
                marker=10,
                title  = mytitle)
        if(N>1)
            for j in 2:(N)
        scatter!((sol[3*N+j,element(i)],sol[4*N+j,element(i)],sol[5*N+j,element(i)]),
                 marker = 5,
                 markerstrokewidth=0,
                 label=planets[j])
            end
        end
    end
   
    gif(anim, filename, fps = 17)

end


#First integrals: energy and angular momentum
export plot_first_integrals
function plot_first_integrals(sol, M::Array{Float64,1}, filename::String, planets::Array{String,1})
    #H(x, y, z, vx, vy, vz, M) need to build the x,y,z and vx,vy,vz arrays
    #initial_first_integrals[1]
    vx=sol[1:N,:]
    vy=sol[N+1:2*N,:]
    vz=sol[2*N+1:3*N,:]
    x= sol[3*N+1:4*N,:]
    y= sol[4*N+1:5*N,:]
    z= sol[5*N+1:6*N,:]
    plot(sol.t, H(x, y, z, vx, vy, vz, M)[1] .- H(x, y, z, vx, vy, vz, M)',
         lab="Energy variation",
         title="First Integrals",
         xlabel="t")
   plot!(sol.t, L(x, y, z, vx, vy, vz, M)[1].-L(x, y, z, vx, vy, vz, M)', lab="Angular momentum variation")

    savefig(filename)
   # return    pot(x, y, z, M)
end



end # module
