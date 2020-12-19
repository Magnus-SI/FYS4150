import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams.update({'font.size': 14})

def comp_plots():
    """
    Compares analytical solution with numerical solutions as function of x
    at certain times.
    """

    dxvals = [0.1, 0.01]
    plotinfo = {'Efor': {'color': 'blue', 'lstyle': ['solid', 'dotted']},
                'imp1': {'color': 'green', 'lstyle': ['solid', 'dotted']},
                'imp2': {'color': 'red', 'lstyle': ['solid', 'dotted']},
                'exct': {'color': 'black', 'lstyle': ['solid', 'dotted']}}

    for dx in dxvals:
        xvals = np.arange(0, 1+1e-15, dx)
        dat = pd.read_csv("data/mcdep_dx_%i.csv"%np.log10(dx),
                          index_col = 0, header = None).T
        plt.figure()
        t1 = float(dat.keys()[0].split("_")[-1])
        for head in dat.keys():
            label = head.split("_")
            t = float(label[-1])
            pinf = plotinfo[label[0][-4:]]
            u = dat[head].values
            if label[0][-4:] == 'exct':
                plt.plot(xvals, u, label = 't = %.1f'%t, color = pinf['color'],
                         linestyle = pinf['lstyle'][t==t1])
            else:
                plt.plot(xvals, u, color = pinf['color'],
                         linestyle = pinf['lstyle'][t==t1])
            plt.title(r"$\Delta x = %.2f, \quad \Delta t = 0.4 \Delta x^2$"%dx)
        plt.plot(xvals, xvals, linestyle = 'dashed', label = r"$t = \infty$")

        plt.legend()
        plt.xlabel("x")
        plt.ylabel("u(x, t_i)")
        plt.show()

        plt.figure()
        plt.title(r"$\Delta x = %.2f, \quad \Delta t = 0.4 \Delta x^2$"%dx)
        lstyles = ['solid', 'dotted']
        for i in range(2):
            u_EF = dat[dat.keys()[4*i]].values
            u_EB = dat[dat.keys()[4*i+1]].values
            u_CN = dat[dat.keys()[4*i+2]].values
            u_exact = dat[dat.keys()[4*i+3]].values

            plt.plot(xvals, u_EF - u_exact, label = r'$\epsilon_{\rm{EF}}, t_%i$'%(i+1),
            color = 'blue', linestyle = lstyles[i])
            plt.plot(xvals, u_EB - u_exact, label = r'$\epsilon_{\rm{EB}}, t_%i$'%(i+1),
            color = 'green', linestyle = lstyles[i])
            plt.plot(xvals, u_CN - u_exact, label = r'$\epsilon_{\rm{CN}}, t_%i$'%(i+1),
            color = 'red', linestyle = lstyles[i])
        plt.xlabel("x")
        plt.ylabel(r"Difference to analytic")
        plt.legend()
        plt.show()

def dt_errplots():
    """
    Plots the error as it depends on dx and dt for different end times and
    the different numerical schemes.
    """
    dxvals = [0.1, 0.01]
    tvals = [0.1, 0.3]
    plotinfo = {'Efor': {'color': 'blue', 'label': "EF, "},
                'imp1': {'color': 'green', 'label': "EB, "},
                'imp2': {'color': 'red', 'label': "CN, "}}
    lstyles = ['solid', 'dashed']

    for t_end in tvals:
        plt.figure()
        plt.title("Global error at t = %.1f"%(t_end))
        for dx, lstyle in zip(dxvals, lstyles):
            dat = pd.read_csv("data/dterr_dxdt_%.3f_%.3f.csv"%(np.log10(dx), t_end))
            dtm = dat['dtm'].values
            for label in dat.keys()[1:]:
                abs_err = dat[label]
                pinf = plotinfo[label]
                plt.plot(np.log2(dtm), np.log10(abs_err), label = pinf['label'] + r"$\Delta x = %.2f$"%dx,
                color = pinf['color'], linestyle = lstyle)
        plt.plot(np.log2(dtm),  np.log10(1e-3 *dtm), color = 'orange', linestyle = 'dotted')
        #plt.plot(np.log10(dtm),  np.log10(1e-4 *dtm**2), color = 'orange', linestyle = 'dotted')
        plt.axvline(-1, color = 'black', linestyle = 'dotted')
        plt.legend()
        plt.xlabel("log2($\Delta t$ multiplier)")
        plt.ylabel("$log(\epsilon_{mean})$")

        plt.ylim(-7, -1)
        plt.show()

def plot2D():
    """
    Plots the 2d case with boundaries = 1 in the upper right corner.
    """
    dx = 0.01
    dat = pd.read_csv("data/diff2D_dx_%.2f.csv"%np.log10(dx),
                      index_col = 0, header = None).T
    for head in dat.keys():
        fig1, ax1 = plt.subplots()
        vals = dat[head].values
        nx = int(np.sqrt(len(vals)))
        X,Y = np.meshgrid(np.linspace(0,1, nx), np.linspace(0,1, nx))
        u = vals.reshape(nx, nx)
        ax = ax1.contourf(X, Y, u, 50, cmap = 'coolwarm')
        fig1.colorbar(ax)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(r"$\Delta x = %.2f, \quad \Delta t = 0.4 \Delta x^2, \quad$ timesteps = %i"%(dx, head))

def err2D():
    """
    Calculates and prints 2d error shown in the table in the report.
    """
    t_end = 1.0
    dat = pd.read_csv("data/err2D_tend_%.2f.csv"%(t_end))
    dat2 = pd.read_csv("data/err2D_local.csv")
    dx = dat['dx'].values
    dt = dat['dt'].values
    errs = dat['err'].values
    dx2 = dat2['dx'].values
    dt2 = dat2['dt'].values
    errs2 = dat2['err'].values
    print("Global errors:")
    for i in range(len(errs)):
        print("dx = %.3f, dt = %.3e, err = %.4e"%(dx[i], dt[i], errs[i]))
    print("\nLocal errors:")
    for i in range(len(errs2)):
        print("dx = %.3f, dt = %.3e, err = %.4e"%(dx2[i], dt2[i], errs2[i]))

def plotlitosphere():
    """
    Plots and saves figures for all litosphere data included in the report
    + some additional figures.
    """
    dx = 0.025
    nx = int(1.5/dx + 1)
    ny = int(1/dx + 1)
    X,Y = np.meshgrid(np.linspace(0,1.5*120, nx), np.linspace(0,1*120, ny))
    #print(nx, ny)
    initdat = pd.read_csv("data/iso2D_stab.csv", index_col = 0, header = None).T
    bvals = initdat[initdat.keys()[0]].values.reshape(ny,nx)

    fig1, ax1 = plt.subplots()
    plt.title("Equilibrium temperature")
    ax1.set_xlabel("x [km]")
    ax1.set_ylabel("depth [km]")
    cont = ax1.contourf(X.T, Y.T, bvals.reshape(ny, nx).T * 1200, 20, cmap = 'plasma')
    cbar = fig1.colorbar(cont)
    cbar.ax.set_ylabel(r"Temperature $[\rm{}\circ C}]$")
    plt.savefig("figures/equib.pdf")

    dat = pd.read_csv("data/iso2D_dx_%.2f.csv"%np.log10(dx),
                      index_col = 0, header = None).T
    #print(dat.keys())
    #bvals = dat[dat.keys()[0]].values.reshape(ny,nx)
    val1 = dat[dat.keys()[-2]].values.reshape(ny,nx)
    val2 = dat[dat.keys()[-1]].values.reshape(ny,nx)
    #plt.figure()
    #plt.imshow(val1 - bvals)
    #plt.figure()
    #plt.imshow(val2 - bvals)

    for i,head in enumerate(dat.keys()):
        fig1, ax1 = plt.subplots()
        plotinfo = head.split("_")
        case = int(plotinfo[0][-1])
        t_end = float(plotinfo[-1][0:])

        plt.title("case %i, t = %.2f Gyr"%(case, t_end/0.01877))
        vals = dat[head].values
        #bvals = dat[dat.keys()[3* (i//3)]].values.reshape(ny,nx)
        u = vals.reshape(ny, nx) - bvals
        #plt.imshow(u)
        ax = ax1.contourf(X.T,Y.T,u.T * 1200, 20, cmap = 'plasma')
        ax1.set_xlabel("x [km]")
        ax1.set_ylabel("depth [km]")
        cbar = fig1.colorbar(ax)
        cbar.ax.set_ylabel(r"$\Delta T$ $[\rm{^\circ C}]$")
        plt.savefig("figures/c%it%.5f.pdf"%(case, t_end))

def plotQ():
    """
    Plots the value of Q which can be saved instead or in addition to T
    in isomain by using the write_Q function in diffusion.cpp.
    Not used in the report, only used for bug-testing.
    """
    dx = 0.025
    nx = int(1.5/dx + 1)
    ny = int(1/dx + 1)
    dat = pd.read_csv("data/Q2D_dx_%.2f.csv"%np.log10(dx),
                      index_col = 0, header = None).T
    X,Y = np.meshgrid(np.linspace(0,1, nx), np.linspace(0,1.5, ny))
    for head in dat.keys():
        plt.figure()
        plt.title(head)
        vals = dat[head].values# - basevals
        u = vals.reshape(ny, nx)
        #plt.imshow(u)
        print(X.shape, Y.shape, u.shape)
        plt.contourf(X, Y, u)

if __name__ == "__main__":
    #Run and plot all results in the report + a few more.
    #Don't run all at the same time, too many plots:

    comp_plots()
    dt_errplots()
    plot2D()
    err2D()
    #plotlitosphere()
