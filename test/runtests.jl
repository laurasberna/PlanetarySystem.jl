using PlanetarySystem
using RecursiveArrayTools
using Test

# The package is tested on random Solar-like systems with a central start initialized at x=(0,0,0) and v=(0,0,0).
# The test verifies the conservation of energy and angular momentum along the evolution, within chosen tolerance.

tol=10^-7 # tolerance

@time @testset "Energy and angular momentum conservation" begin

    for j in 1:10  # loop over different random systems
        
        # initial conditions:
        N=5 # number of bodies
        M=rand(0.:0.000001:0.0001,N) # mass
        M[1]=1. # central star mass
        pos_x=rand(-30.:0.0001:30.,N) # initial positions
        pos_y=rand(-30.:0.0001:30.,N)
        pos_z=rand(-30.:0.0001:30.,N)
        pos_x[1]=0. # central star 
        pos_y[1]=0.
        pos_z[1]=0.
        pos = ArrayPartition(pos_x,pos_y,pos_z)
        vel_x=rand(-0.01:0.00001:0.01,N) # initial velocities
        vel_y=rand(-0.01:0.00001:0.01,N)
        vel_z=rand(-0.01:0.00001:0.01,N)
        vel = ArrayPartition(vel_x,vel_y,vel_z)
        vel_x[1]=0. # central star
        vel_y[1]=0.
        vel_z[1]=0.
        tspan = (0.,200_000.) # time interval
        
        # solution:
        sol=NBsolution(M, vel, pos, tspan)
        # solution partition:
        vx=sol[1:N,:]
        vy=sol[N+1:2*N,:]
        vz=sol[2*N+1:3*N,:]
        x= sol[3*N+1:4*N,:]
        y= sol[4*N+1:5*N,:]
        z= sol[5*N+1:6*N,:]

        # conservation test:
        if(length(sol.t)>=5000) # test step check
            
            for i in 1:5000:length(sol.t) # loop over the time steps; test step=500
                @test (H(x, y, z, vx, vy, vz, M)[1] -  H(x, y, z, vx, vy, vz, M)[i])/H(x, y, z, vx, vy, vz, M)[1] ≈ 0 atol=tol #energy
                @test (Lx(x, y, z, vx, vy, vz, M)[1] - Lx(x, y, z, vx, vy, vz, M)[i])/Lx(x, y, z, vx, vy, vz, M)[1] ≈ 0 atol=tol #angular momentum x
                @test (Ly(x, y, z, vx, vy, vz, M)[1] - Ly(x, y, z, vx, vy, vz, M)[i])/Ly(x, y, z, vx, vy, vz, M)[1]  ≈ 0 atol=tol #angular momentum y
                @test (Lz(x, y, z, vx, vy, vz, M)[1] - Lz(x, y, z, vx, vy, vz, M)[i])/Lz(x, y, z, vx, vy, vz, M)[1] ≈ 0 atol=tol #angular momentum z
            end
            
        else
            println("reduce test step!") # alert
        end

     end
    
end


