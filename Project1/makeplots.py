import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def plot_sol(method):
    """
    If method is lin:
    Plots the values of the approximate solution v against the analytical solution
    u as a function of x for values of n = 10, 100, 1000.
    Saves the result in figure lin.pdf.

    If method is log:
    plots the relative error between u and v in a logarithmic plot as a function of x.
    Saves the result in figure log.pdf.

    The used values are stored in valuesn.csv for the values of n.
    """
    plt.figure()
    for N in ('1', '2', '3'):
        df = pd.read_csv('vals/values' + N + '.csv')
        x = df['x'].values
        v = df['v'].values
        if method == 'lin':
            #Plotting solution for different N
            plt.plot(x, v, label = r'n = 10$^{%s}$'%N)
            if N == '3':
                u = df['u'].values
                plt.plot(x, u, label = 'closed form sol.')
                plt.xlabel(r"$x$")
                plt.ylabel(r"$v$")
        elif method == 'log':
            #Plotting relative error logarithmically
            u = df['u'].values
            plt.semilogy(x, np.abs((u-v)/u), label = 'N = %s'%N)
    plt.legend()
    plt.savefig(method+".pdf")

def plot_maxeps():
    """
    Plots the maximum value of log10 of the relative errors from the results
    from the optimized method (arbitrary choice as v is the same for all methods,
    but optimized is the fastest way to generate these results).
    The values used are stored in epsclock0.csv.
    The resulting figure is saved in eps.pdf
    """
    plt.figure()
    data = pd.read_csv('vals/epsclock0.csv')
    N = data['expv'].values
    eps_max = data['maxeps'].values
    plt.plot(N, eps_max, label='Measured error')
    plt.xlabel(r"log$_{10}$(n)")
    plt.ylabel(r"log$_{10}$($\epsilon_{max}$)")
    plt.plot(N, -2*N+eps_max[0]+2, ls='--', label=r'$\mathcal{O}(n^{-2})$')
    plt.legend()
    plt.savefig("eps.pdf")

def plot_time():
    """
    Compares the computational time of the different algorithms and plots
    them as a function of grid points n.
    The resulting figure is saved as time.pdf
    """
    plt.figure()
    optimized = pd.read_csv('vals/epsclock0.csv')
    N_o = optimized['expv'].values
    time_o = optimized['time'].values

    general = pd.read_csv('vals/epsclock1.csv')
    N_g = general['expv'].values
    time_g = general['time']
    time_LU = general['LUtime']

    plt.plot(N_o, np.log10(time_o), label='Optimized algorithm')
    plt.plot(N_g, np.log10(time_g), label='General algorithm')
    plt.plot(N_g, np.log10(time_LU), label='LU decomposition')
    plt.plot(N_o, N_o+np.log10(time_o)[0]-2, linestyle = '--', label=r'$\mathcal{O}(n)$')
    plt.plot(N_g[:-3], 3*N_g[:-3]+np.log10(time_o)[0]-4, linestyle = '--', label=r'$\mathcal{O}(n^3)$')
    plt.xlabel(r"log$_{10}$(n)")
    plt.ylabel(r"log$_{10}$(runtime)")
    plt.legend()
    plt.savefig("time.pdf")

if __name__ == "__main__":
    "Generate the plots used in the report:"
    plot_sol("lin")
    plot_maxeps()
    plot_time()
