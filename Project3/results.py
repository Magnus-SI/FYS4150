import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import constants, units
from mpl_toolkits.mplot3d import Axes3D

def mdimarr(data, nt, n):
    x = data["x"]*units.m.to("au")
    y = data["y"]*units.m.to("au")
    z = data["z"]*units.m.to("au")
    vx = data["vx"]
    vy = data["vy"]
    vz = data["vz"]

    r = np.zeros((3, nt, n))
    v = np.copy(r)
    for i in range(n):
        r[0, :, i] = x[i::n]
        r[1, :, i] = y[i::n]
        r[2, :, i] = z[i::n]
        v[0, :, i] = vx[i::n]
        v[1, :, i] = vy[i::n]
        v[2, :, i] = vz[i::n]

    return r, v


def betaplots():
    nt = int(1e4)
    n = 2
    betas = np.array([2, 2.001, 2.01, 2.1, 3])
    plt.figure()
    for beta in betas:
        data = pd.read_csv("data/earth_sun2_%.2f_4.00.txt"%beta)
        r,v = mdimarr(data, nt, n)
        plt.plot(r[0,:,1], r[1,:,1], label = r"$\beta = %.2f$"%beta)
    plt.legend()
    plt.show()

def escvelplots():
    pass

def jupiterplots():
    pass

def solarplots():
    pass

def mercury_recession():
    """
    Plots mercury orbit and calculates precession?
    """
    data = pd.read_csv("data/mercury_2_2.00_5.00.txt")
    nt = int(1e5)
    n = 2
    r, v = mdimarr(data, nt, n)
    pass

if __name__ == "__main__":
    #call functions here as wanted
    pass
