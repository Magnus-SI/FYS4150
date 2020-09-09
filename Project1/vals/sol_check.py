"""
Script for comparing solutions from
general, optimized, and LU algorithms
for N = 10, 100, 1000
"""

import pandas as pd
import numpy as np

for exponent in [1, 2, 3]:
    sol = pd.read_csv('values'+str(exponent)+'.csv')
    gen = sol['v'].values
    LU = sol['LUv'].values
    if np.any(np.abs((gen - LU)) > 0):
        print("General algorithm and LU algorithm don't match")
        print("Maximum difference of",np.max(np.abs(gen-LU)),"in n = 10^%i" %exponent)
    else:
        print("LU and general algorithms yielded exact same result at n = 10^%i" %exponent)
