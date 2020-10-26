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
    nt = int(4e4)
    n = 2
    betas = np.array([2, 2.001, 2.01, 2.1, 2.2])
    plt.figure()
    for beta in betas:
        data = pd.read_csv("data/elliptical_earth_sun2_%.2f_4.60.txt"%beta)
        r,v = mdimarr(data, nt, n)
        plt.plot(r[0,:,1], r[1,:,1], label = r"$\beta = %.2f$"%beta)
    plt.legend()
    plt.show()
    plt.figure()
    for beta in betas:
        data = pd.read_csv("data/earth_sun2_%.2f_4.60.txt"%beta)
        r,v = mdimarr(data, nt, n)
        plt.plot(r[0,:,1], r[1,:,1], label = r"$\beta = %.2f$"%beta)
    plt.legend()
    plt.show()

def escvelplots():
    pass

def jupiterplots():
    nt = 20000
    n = 3
    massms = np.array([1, 1000])

    for i, title in enumerate(["Fixed sun", "Non-fixed sun"]):
        plt.figure(i)
        plt.axis([-10, 10, -10, 10])
        plt.gca().set_aspect("equal", adjustable = 'box')
        plt.title(title)

    for m, lstyle in zip(massms, ["-", "dotted"]):
        fdata =  pd.read_csv("data/m%.1fjupiter3_2.00_4.30.txt"%m)
        nfdata = pd.read_csv("data/nfm%.1fjupiter3_2.00_4.30.txt"%m)
        r1,v1 = mdimarr(fdata, nt, n)
        r2,v2 = mdimarr(nfdata, nt, n)
        for i, pname, c in zip(range(0,3), ["sun", "earth", "jupiter"], ["red", "blue", "orange"]):
            plt.figure(0)
            plt.plot(r1[0, :, i], r1[1, :, i], linestyle = lstyle,
                     label = "%s, multiplier %i"%(pname,m), color = c)
            plt.figure(1)
            plt.plot(r2[0, :, i], r2[1, :, i], linestyle = lstyle,
                     label = "%s, multiplier %i"%(pname,m), color = c)

        plt.legend()
        plt.figure(0)
        plt.legend()

def solarplots():
    nt = 20000
    n = 10
    data = pd.read_csv("data/solar%i_2.00_%.2f.txt"%(n, np.log10(nt)))
    r,v = mdimarr(data, nt, n)
    plt.figure()
    for i in range(n):
        plt.plot(r[0, :, i], r[1, :, i], label = "planet %i"%i)
    plt.legend()
    plt.show()

    energy_data = pd.read_csv("tot_en.csv")
    t = energy_data["      timestep"]
    PE = energy_data["PE"]
    KE = energy_data["KE"]
    plt.figure()
    plt.plot(t, PE, label = "PE")
    plt.plot(t, KE, label = "KE")
    plt.yscale("log")
    plt.legend()
    plt.show()


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
