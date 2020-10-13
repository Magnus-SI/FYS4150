import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

data = pd.read_csv("solar.txt")

x = data["x"]
y = data["y"]
z = data["z"]
vx = data["vx"]
vy = data["vy"]
vz = data["vz"]


"""fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x,y,z)
plt.show()"""

plt.plot(x,y,'.')
plt.show()
