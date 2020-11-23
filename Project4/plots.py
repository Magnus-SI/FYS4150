import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.stats import linregress
plt.rcParams.update({'font.size': 14})

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
    Tmin = 2.0
    Tmax = 2.3
    Tvals = np.linspace(Tmin, Tmax, Tcount)
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
    Tcount = 16
    Lcount = 4
    Tmin = 2.0
    Tmax = 2.3
    Ls = np.array([40, 60, 80, 100])
    Tvals = np.linspace(Tmin, Tmax, Tcount)
    Cv = np.zeros((Lcount, Tcount))
    for i,T in enumerate(Tvals):
        data = pd.read_csv("data/T%.6fmultiL.csv"%T)
        Cv[:, i] = data['Cv'].values

    Tnew = np.linspace(Tmin, Tmax, 1501)
    Cvfunc = interp1d(Tvals, Cv, axis = 1, kind ='cubic')
    Cvnew = Cvfunc(Tnew)
    Tcs = Tnew[np.argmax(Cvnew, axis = 1)]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("Cubic spline interp of $C_v(T)$")
    for i in range(len(Ls)):
        line, = ax.plot(Tvals, Cv[i], 'x', label = "%ix%i"%(Ls[i], Ls[i]))
        ax.plot(Tnew, Cvnew[i], linestyle = 'dashed', color = line.get_color())
        ax.axvline(x = Tcs[i], linestyle = 'dotted', color = line.get_color())
    ax.set_xlabel("Temperature")
    ax.set_ylabel(r"$C_v$")
    plt.legend()
    plt.show()



    x = 1/Ls
    x_long = 1/np.linspace(30,110,81)
    slope, intercept, rvalue, pvalue, stderr = linregress(x, Tcs)

    plt.figure()
    plt.plot(x, Tcs, 'rx')
    plt.plot(x_long, intercept + slope*x_long, 'g--',
    label = "$T_C(L=\infty) = %.4f$,  $a = %.1f$,\n$\sigma_{\\rm{err}} = %.2e$"%(intercept, slope, stderr))
    plt.legend()
    plt.title("Linear fit of $T_C(L) = T_C(L=\infty) + aL^{-1}$")
    plt.xlabel(r"$L^{-1}$")
    plt.ylabel("$T_c$")
    plt.show()
    print(stderr)

def timerplots():
    timerdat = pd.read_csv("data/timer.csv")
    Lvals = timerdat['L'].values
    npcs = timerdat['npcs'].values
    t = timerdat['time'].values
    L_count = len(np.unique(Lvals))
    npc_count = len(np.unique(npcs))
    pcrange = np.linspace(1, 8, 71)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("Speedup vs. Predicted $\propto$ 1/(pc_num)\nmcs = 10 000")
    pcL = np.unique(npcs)
    for L, color in zip(np.unique(Lvals), ['blue', 'green']):
        tL = np.zeros(npc_count)
        for j,npc in enumerate(pcL):
            tr_arr = (Lvals == L) * (npcs == npc)
            ts = t[tr_arr]
            tL[j] = np.sum(ts)/len(ts)
        line, = ax.plot(pcL, tL, "x", label = "%ix%i"%(L, L), color = color)
        ax.plot(pcrange, tL[0]/pcrange, linestyle = 'dashed', color = line.get_color())

    plt.legend()
    plt.grid()
    ax.set_xlabel("Processor count")
    ax.set_ylabel("time[s]")
    plt.show()



if __name__ == "__main__":
    #mcplots()
    #TLvar('E', r"$E/\rm{J}$")
    #TLvar('M', r'$\langle | M | \rangle$')
    #TLvar('Cv', r'$C_V / \rm{Jk_B}$')
    #TLvar('chi', r'$\chi$')
    #Tcestim()
    timerplots()
