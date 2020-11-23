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
            for quant in dat.keys()[1:-1]:
                qval = dat[quant].values
                err = np.abs((qval[1:] - qval[:-1])/qval[:-1])
                plt.plot(dat['mc'].values[1:], err, label  = quant)
            plt.legend()
            plt.yscale("log")
            plt.show()

            #1e:
            plt.figure()
            plt.title("T = %.2f, tol = %.2f"%(T, tol))
            E_equi = dat['E'].values[7000:]
            plt.hist(E_equi, bins = 20)
            plt.show()

            plt.figure()
            plt.title("T = %.2f, tol = %.2f"%(T, tol))
            acptfrac = dat['acptfrac'].values[7000:]
            plt.plot(dat['mc'].values, dat['acptfrac'].values)
            plt.show()

def TLvar(quantlabel, plotlabel):
    Tcount = 16
    Lcount = 4
    Tmin = 2
    Tmax = 2.3
    Tvals = np.linspace(2, 2.3, Tcount)
    quant = np.zeros((Lcount, Tcount))

    for i,T in enumerate(Tvals):
        data = pd.read_csv("data/T%.6fmultiL.csv"%T)
        quant[:, i] = data[quantlabel].values

    Lvals = data['L'].values

    plt.figure()
    plt.ylabel(plotlabel)
    plt.xlabel(r"$k_B T$")
    for i in range(Lcount):
        plt.plot(Tvals, quant[i,:], 'x', label = "%ix%i"%(Lvals[i], Lvals[i]))
    plt.legend()
    plt.show()

def Tcestim():
    #Gonna use cubic splines here.
    pass

if __name__ == "__main__":
    #mcplots()
    TLvar('E', r"$E/\rm{J}$")
    TLvar('M', r'$\langle | M | \rangle$')
    TLvar('Cv', r'$C_V / \rm{Jk_B}$')
    TLvar('chi', r'$\chi$')
