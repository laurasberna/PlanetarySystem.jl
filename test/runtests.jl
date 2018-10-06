using PlanetarySystem
#using Test
using RecursiveArrayTools
using Plots

plotfilename="mysol.png"
animfilename="mysol.gif"
fifilename="firstintegrals.png"
mytitle="Planetary System (outer bodies of Solar system)"

#System parameters:
#masses and names
M = [1.00000597682, 0.000954786104043, 0.000285583733151, 0.0000437273164546, 0.0000517759138449]#, 1/1.3e8]
invM = inv.(M)
planets = ["Sun", "Jupiter", "Saturn", "Uranus", "Neptune"] #, "Pluto"
#initial positions
pos_x = [0.0,-3.5023653,9.0755314,8.3101420,11.4707666]#,-15.5387357]
pos_y = [0.0,-3.8169847,-3.0458353,-16.2901086,-25.7294829]#,-25.2225594]
pos_z = [0.0,-1.5507963,-1.6483708,-7.2521278,-10.8169456]#,-3.1902382]
pos = ArrayPartition(pos_x,pos_y,pos_z)
#initial velocities
vel_x = [0.0,0.00565429,0.00168318,0.00354178,0.00288930]#,0.00276725]
vel_y = [0.0,-0.00412490,0.00483525,0.00137102,0.00114527]#,-0.00170702]
vel_z = [0.0,-0.00190589,0.00192462,0.00055029,0.00039677]#,-0.00136504]
vel = ArrayPartition(vel_x,vel_y,vel_z)
#time interval
tspan = (0.,200_000)

#Testing solver
mysol=NBsolution(M, vel, pos, tspan)

#Testing orbit plot
#myplot(mysol,plotfilename,planets,mytitle)

#Testing animation 
#animation(mysol,animfilename,planets,mytitle)

#Initial energy
#@show H(pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,M)

#Plot energy
N=5
#@show length(mysol.t)
#@show length(plot_first_integrals(mysol, M, fifilename, planets))
plot_first_integrals(mysol, M, fifilename, planets)
#@show length(mysol[3*N+1:4*N,:])
#@show mysol[1:N,:]
#x=zeros(N,length(mysol.t))
#for i in 1:N x[i,:]=mysol[3*N+i,:] end
#@show x[1,:]
#=    vx=sol[1:N,1:2]
    vy=sol[N+1:2*N,1:2]
    vz=sol[2*N+1:3*N,1:2]
    x= sol[3*N+1:4*N,1:2]
    y= sol[4*N+1:5*N,1:2]
z= sol[5*N+1:6*N,1:2]
@show H(x, y, z, vx, vy, vz, M)
@show H(x, y, z, vx, vy, vz, M)[1]=#
#@show vx.*vx
#@show sum(M.*vx.*vx,dims=1)


#=Testing write to file  (not working)
filename="sol.h5"
rm(filename, force=true)
writefile(mysol,filename)

#Testing read from file
mysol2=readfile(filename)
#@test isequal(mysol2,mysol)=#
