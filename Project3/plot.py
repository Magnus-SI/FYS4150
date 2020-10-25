import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import constants, units
from mpl_toolkits.mplot3d import Axes3D

nt = 1000
n = 10

data = pd.read_csv("solar10_2.00_3.00.txt")

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
#plt.ylim([-2,2])
#plt.xlim([-2, 2])
plt.axis("equal")
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for j in range(n):
    ax.plot(X[0::n,j],Y[0::n,j],Z[0::n,j])

plt.show()

mercury = pd.read_csv("mercury.txt")

x_M = mercury["x"]*units.m.to("au")
y_M = mercury["y"]*units.m.to("au")
z_M = mercury["z"]*units.m.to("au")
vx_M = mercury["vx"]
vy_M = mercury["vy"]
vz_M = mercury["vz"]

X_M = np.zeros((nt, 2))
Y_M = X_M.copy()
Z_M = X_M.copy()

for i in range(2):
    X_M[:,i] = x_M[i::2]
    Y_M[:,i] = y_M[i::2]
    Z_M[:,i] = z_M[i::2]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(X_M[:,1],Y_M[:,1],Z_M[:,1])
plt.show()

plt.plot(X_M[:,0], Y_M[:,0])
plt.plot(X_M[:,1], Y_M[:,1])
plt.axis("equal")
plt.show()
