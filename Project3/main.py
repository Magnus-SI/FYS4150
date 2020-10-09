import numpy as np
import os
import pandas as pd

initial_raw = open("initial_raw.txt", "r")
print("Opening intitial_raw.txt and extracting initial values")

objects = []
mass = []       #10e24 kg
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
        mass.append(line[7:].rstrip())
    elif line[0] == "X":
        X.append(line[3:25])
        Y.append(line[29:51])
        Z.append(line[55:].rstrip())
    elif line[:2] == "VX":
        VX.append(line[3:25])
        VY.append(line[29:51])
        VZ.append(line[55:].rstrip())

initial_raw.close()

N = len(objects)

initial = open("initial.txt","w")
print("Writing initial positions and velocities to file")
for line in range(N):
    initial.write(X[line]+" "+Y[line]+" "+Z[line]+" ")
    initial.write(VX[line]+" "+VY[line]+" "+VZ[line])
    initial.write("\n")
initial.close()

mass_file = open("masses.txt","w")
print("Writing mass for objects to file")
for line in range(N):
    mass_file.write(mass[line]+"\n")
mass_file.close()
