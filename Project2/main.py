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

eig = pd.read_csv("test.csv")
print(eig.keys())

plt.plot(eig["eigenvector1"])
plt.plot(eig["eigenvector2"])
plt.plot(eig["eigenvector3"])
plt.show()
