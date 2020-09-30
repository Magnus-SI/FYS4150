import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
plt.rcParams.update({'font.size': 14})

def plotbeam():
    """
    Plot the buckling beam figure shown in the report.
    """
    plt.figure()
    eigvecs = pd.read_csv('csv_files/beam.csv')
    x = np.linspace(0, 1, 102)
    for i in range(3):
        u = np.zeros(len(x))
        u[1:-1] = eigvecs['eigenvector%i'%(i+1)]
        plt.plot(x, u, label = r'$v_%i$'%(i+1))
    plt.legend()
    plt.show()

def plotelectron2():
    """
    Plot the 2-electron quantum system figure shown in the report.
    """
    plt.figure()
    omega_rs = [0.01, 0.5, 1, 5]
    rhomax = 5
    rhos = np.linspace(0+(rhomax)/101, rhomax-rhomax/101, 100)
    for i in range(4):
        eigvecs = pd.read_csv("csv_files/el2omega%i.csv"%i)
        plt.plot(rhos, eigvecs["eigenvector1"], label = r"$\omega_r = %g$"%omega_rs[i])
    plt.legend()
    plt.xlabel(r"$\rho$")
    plt.show()

def plotiter():
    """
    Plot the computational time and iteration plot shown in the report.
    """
    plt.figure()
    iter = pd.read_csv("csv_files/timer.csv")
    print(iter.keys)
    log10N = iter["        log10N"]
    iters = iter["iters"]
    Jtime = iter["Jtime"]
    plt.plot(log10N, np.log10(iters), color='b', label='y = iterations')
    #plt.plot(log10N, np.log10((10**log10N)**2 - 10**log10N), ls="--")
    #plt.plot(log10N, np.log10((10**log10N)**2), ls="--")
    plt.plot(log10N, np.log10(1.4*((10**log10N)**2 - 10**log10N)), ls="--", color='b',\
        label=r"$y = N^2 - N$")
    plt.xlabel(r"$log_{10}(N)$")
    plt.ylabel(r"$\log_{10}$(y)")
    plt.plot(log10N, np.log10(Jtime*1000), color='g', label='y = runtime [ms]')
    plt.plot(log10N, np.log10((10**log10N)**4 - (10**log10N)**3)-4.5 , color='g', ls='--',\
        label=r"$y = N^4 - N^3$")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    plotbeam()
    plotelectron2()
    plotiter()
