import numpy as np
import os
import pandas as pd
from astropy import constants, units
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

initial_raw = open("initial_raw.txt", "r")
print("Opening intitial_raw.txt and extracting initial values")

n = int(input("How many planets to include (Pluto is a plantet!): ")) + 1
nt = int(input("How many timesteps: "))

objects = []
mass = []       #1e24 kg
X = []          #AU
Y = []
Z = []
VX = []         #AU/day
VY = []
VZ = []

for line in initial_raw:
    if line[0] == "#":
        objects.append(line[1:].rstrip())
    elif line[:4] == "Mass":
        mass.append(float(line[7:].rstrip()))
    elif line[0] == "X":
        X.append(float(line[3:25]))
        Y.append(float(line[29:51]))
        Z.append(float(line[55:].rstrip()))
    elif line[:2] == "VX":
        VX.append(float(line[3:25]))
        VY.append(float(line[29:51]))
        VZ.append(float(line[55:].rstrip()))

initial_raw.close()

N = len(objects)

#converting to SI units
mass = np.array(mass)*1e24
X = np.array(X)*constants.au.to("m").value
Y = np.array(Y)*constants.au.to("m").value
Z = np.array(Z)*constants.au.to("m").value
VX = np.array(VX)*constants.au.to("m").value/86400
VY = np.array(VY)*constants.au.to("m").value/86400
VZ = np.array(VZ)*constants.au.to("m").value/86400

initial = open("initial.txt","w")
print("Writing initial positions and velocities to file")
for line in range(n):
    initial.write(str(X[line])+" "+str(Y[line])+" "+str(Z[line])+" ")
    initial.write(str(VX[line])+" "+str(VY[line])+" "+str(VZ[line]))
    initial.write("\n")
initial.close()

mass_file = open("masses.txt","w")
print("Writing mass for objects to file")
for line in range(n):
    mass_file.write(str(mass[line])+"\n")
mass_file.close()

main = "main.cpp"
super = "solar_system.cpp"
exe = "main.out"
os.system("echo compiling programs...")
compile = " ".join(["c++", "-o", exe, main, super])
os.system(compile)
os.system("./"+exe+" "+str(n)+" "+str(nt))

data = pd.read_csv("solar.txt")

x = data["x"]*units.m.to("au")
y = data["y"]*units.m.to("au")
z = data["z"]*units.m.to("au")
vx = data["vx"]
vy = data["vy"]
vz = data["vz"]

X = np.zeros((nt,n))
Y = np.copy(X)
Z = np.copy(Y)

for i in range(n):
    X[:,i] = x[i::n]
    Y[:,i] = y[i::n]
    Z[:,i] = z[i::n]

for k in range(n):
    plt.plot(X[:,k], Y[:,k])
plt.axis("equal")
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for j in [0,1,3,4,6]:
    ax.plot(X[0::10,j],Y[0::10,j],Z[0::10,j])    
plt.show()

"""plt.plot(x,y,'.')
plt.axis("equal")
plt.show()"""
