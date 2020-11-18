import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def mcplots():
    temps = [1.0, 2.4]
    tols = [0.0, 0.5]
    for T in temps:
        for tol in tols:
            dat = pd.read_csv("data/mcdep_20_%.2f_%.2f.csv"%(T, tol))
            plt.figure()
            plt.title("T = %.2f, tol = %.2f"%(T, tol))
            for quant in dat.keys()[1:]:
                qval = dat[quant].values
                err = np.abs((qval[1:] - qval[:-1])/qval[:-1])
                plt.plot(dat['mc'].values[1:], err, label  = quant)
            plt.legend()
            plt.yscale("log")
            plt.show()

if __name__ == "__main__":
    mcplots()
