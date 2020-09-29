import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

"""
if not os.path.exists("main.out"):
    os.system("c++	-o	main.out	main.cpp Jacobi_rotation.cpp	-larmadillo")
if not os.path.exists("test-main.out"):
    os.system("c++ -o tests.out tests-main.cpp test-functions.cpp Jacobi_rotation.cpp -larmadillo")
"""
#
# eig = pd.read_csv("test.csv")
# print(eig.keys())
#
# plt.plot(eig["eigenvector1"])
# plt.plot(eig["eigenvector2"])
# plt.plot(eig["eigenvector3"])
# plt.show()

def plotelectron2():
    plt.figure()
    omega_rs = [0.01, 0.5, 1, 5]
    rhomax = 4
    rhos = np.linspace(0+(rhomax)/101, rhomax-rhomax/101, 100)
    for i in range(4):
        eigvecs = pd.read_csv("el2omega%i.csv"%i)
        plt.plot(rhos, eigvecs["eigenvector1"], label = r"$\omega_r = %g$"%omega_rs[i])
    plt.legend()
    plt.xlabel(r"$\rho$")
    plt.show()
