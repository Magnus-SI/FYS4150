import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv("solar.txt")

x = data["x"]
y = data["y"]
z = data["z"]
vx = data["vx"]
vy = data["vy"]
vz = data["vz"]

plt.plot(x,y,'.')
plt.axis("equal")
plt.show()
