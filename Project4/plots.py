import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def mcplots():
    dat = pd.read_csv("f1.csv")
    plt.figure()
    expvals = [0, 8, 256, 16]
    for quant, expval in zip(dat.keys()[1:], expvals):

        err = np.abs(dat[quant] - expval)
        plt.plot(dat['mc'], err, label  = quant)
    plt.legend()
    plt.yscale("log")
    plt.show()

if __name__ == "__main__":
    mcplots()
