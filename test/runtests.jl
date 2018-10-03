using PlanetarySystem
#using Test
using RecursiveArrayTools
using Plots

plotfilename="mysol.png"
animfilename="mysol.gif"
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

#Testing plot
#myplot(mysol,planets,plotfilename,mytitle)

#=Testing write to file
filename="sol.h5"
rm(filename, force=true)
writefile(mysol,filename)

#Testing read from file
mysol2=readfile(filename)
#@test isequal(mysol2,mysol)=#

#Testing animation
#animation(mysol,animfilename)

(mysol[16,1],mysol[21,1],mysol[26,1])
