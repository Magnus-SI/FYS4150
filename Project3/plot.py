import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import constants, units
from mpl_toolkits.mplot3d import Axes3D

nt = 1000
n = 2

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
#plt.ylim([-2,2])
#plt.xlim([-2, 2])
plt.axis("equal")
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for j in range(n):
    ax.plot(X[0::n,j],Y[0::n,j],Z[0::n,j])

plt.show()
