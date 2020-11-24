import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.stats import linregress
plt.rcParams.update({'font.size': 14})

def TLvar(quantlabel, plotlabel, title):
    """
    Plot quantity for different lattice sizes as function of temperature
    """
    Tcount = 18
    Lcount = 4
    Tmin = 2.0
    Tmax = 2.34
    Tvals = np.linspace(Tmin, Tmax, Tcount)
    quant = np.zeros((Lcount, Tcount))

    for i,T in enumerate(Tvals):
        data = pd.read_csv("data/T%.6fmultiL.csv"%T)
        quant[:, i] = data[quantlabel].values

    Lvals = data['L'].values

    plt.figure(figsize = [7,5.3])
    plt.title(title)
    plt.ylabel(plotlabel)
    plt.xlabel(r"$k_B T$")
    for i in range(Lcount):
        plt.plot(Tvals, quant[i,:], 'x', label = "%ix%i"%(Lvals[i], Lvals[i]))
    plt.legend()
    plt.savefig("figs/diff_" + quantlabel + "s.pdf")
    plt.show()

def Tcestim():
    """
    Estimates the crticial temperature
    """
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
    slope, intercept, rval, pval, stderr = linregress(x, Tcs)

    T_hat = intercept + slope*x
    n = len(x)
    x_mean = np.mean(x)
    SSxx = np.sum((x - x_mean)**2)
    SSE = np.sum((Tcs - T_hat)**2)
    syx2 = SSE/(n-2)
    std_b = np.sqrt(syx2 * (1/n + x_mean**2/SSxx))


    plt.figure()
    plt.plot(x, Tcs, 'rx')
    plt.plot(x_long, intercept + slope*x_long, 'g--',
    label = "$T_C(L=\infty) = %.4f$,  $a = %.1f$,\n$\sigma_{\\rm{b}} = %.2e$"%(intercept, slope, std_b))
    plt.legend()
    plt.title("Linear fit of $T_C(L) = T_C(L=\infty) + aL^{-1}$")
    plt.xlabel(r"$L^{-1}$")
    plt.ylabel("$T_c$")
    plt.show()



def timerplots():
    """
    Plots data for timing of the parallelization.
    """
    timerdat = pd.read_csv("data/timer.csv")
    Lvals = timerdat['L'].values
    npcs = timerdat['npcs'].values
    t = timerdat['time'].values
    L_count = len(np.unique(Lvals))
    npc_count = len(np.unique(npcs))
    pcrange = np.linspace(1, 8, 71)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("Speedup vs. Predicted $\propto$ 1/#cores\nmcs = 10 000")
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
    ax.set_xlabel("Core count")
    ax.set_ylabel("time[s]")
    plt.show()



if __name__ == "__main__":
    #mcplots()
    TLvar('E', r"$\langle E/N \rangle /\rm{J}$", "Mean energy per spin")
    TLvar('M', r'$\langle | M/N | \rangle$', "Mean absolute magnetisation per spin")
    TLvar('Cv', r'$C_V/N / \rm{Jk_B}$', "Heat capcacity")
    TLvar('chi', r'$\chi/N$', "Magnetic susceptibility")
    Tcestim()
    timerplots()
