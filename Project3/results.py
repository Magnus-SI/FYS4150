import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import constants, units
from mpl_toolkits.mplot3d import Axes3D

plt.rcParams.update({'font.size': 14})


def mdimarr(data, nt, n):
    x = data["x"]*units.m.to("au")
    y = data["y"]*units.m.to("au")
    z = data["z"]*units.m.to("au")
    vx = data["vx"]
    vy = data["vy"]
    vz = data["vz"]

    r = np.zeros((3, nt, n))
    v = np.copy(r)
    for i in range(n):
        r[0, :, i] = x[i::n]
        r[1, :, i] = y[i::n]
        r[2, :, i] = z[i::n]
        v[0, :, i] = vx[i::n]
        v[1, :, i] = vy[i::n]
        v[2, :, i] = vz[i::n]

    return r, v

def err_plots():
    data = pd.read_csv("data/timeacc.csv")
    logN = data["        log10N"]
    eultime = data["eultime"]
    vertime = data["vertime"]
    eulerr = data["eulerr"]
    vererr = data["vererr"]
    plt.figure()
    line1, = plt.plot(logN, np.log10(eultime), label=r"Euler time, $\log (t)$ [s]", color='b')
    line2, = plt.plot(logN, np.log10(vertime), label=r"Verlet time, $\log (t)$ [s]", color='g')
    line3, = plt.plot(logN, eulerr, ls="--", color='b', label="Euler accuracy")
    line4, = plt.plot(logN, vererr, ls='--', color='g', label="Verlet accuracy")
    line5, = plt.plot(logN, logN-7, ls='dashdot', color='r', label=r"$\mathcal{O}(N)$")
    plt.xlabel(r"$\log_{10}(N)$")
    # Create a legend for the first line.
    first_legend = plt.legend(handles=[line2, line1, line5], loc='lower right')

    # Add the legend manually to the current Axes.
    ax = plt.gca().add_artist(first_legend)

    # Create another legend for the second line.
    plt.legend(handles=[line4, line3], loc='center left')
    plt.xticks([3,4,5,6])
    plt.tight_layout()
    plt.savefig("figures/err_plots.pdf")

def betaplots():
    nt = int(1e5)
    n = 2
    betas = np.array([2, 2.01, 2.1, 2.2, 2.5])
    plt.figure()
    for i,beta in enumerate(betas):
        data = pd.read_csv("data/elliptical_earth_sun2_%.3f_%.3f.txt"%(beta, np.log10(nt)))
        r,v = mdimarr(data, nt, n)
        if i == 0:
            plt.plot(r[0,:int(((i+1)/6)**4*nt),1][::int(((i+1)/6)**2*100)], r[1,:int(((i+1)/6)**4*nt),1][::int(((i+1)/6)**2*100)], label = r"$\beta = %.3f$"%beta\
                , color ='k', linewidth=4)
        else:
            plt.plot(r[0,:int(((i+1)/6)**4*nt),1][::int(((i+1)/6)**2*100)], r[1,:int(((i+1)/6)**4*nt),1][::int(((i+1)/6)**2*100)], label = r"$\beta = %.3f$"%beta\
                , alpha = 0.7)
    plt.legend()
    plt.axis("equal")
    plt.xlabel(r"$x$ [AU]")
    plt.ylabel(r"$y$ [AU]")
    plt.tight_layout()
    plt.savefig("figures/beta_elliptical.pdf")
    plt.figure()
    betas = np.array([2, 2.1, 2.5, 3])
    for i,beta in enumerate(betas):
        data = pd.read_csv("data/earth_sun2_%.3f_%.3f.txt"%(beta, np.log10(nt)))
        r,v = mdimarr(data, nt, n)
        plt.plot(r[0,:,1][::int(((i+1)/6)**2*100)], r[1,:,1][::int(((i+1)/6)**2*100)] + 0.1*i, label = r"$\beta = %.3f$"%beta)
    plt.axis("equal")
    plt.xlabel(r"$x$ [AU]")
    plt.ylabel(r"$y$ [AU]")
    plt.legend()
    plt.tight_layout()
    plt.savefig("figures/beta_circular.pdf")


def escvelplots():
    nt = int(4e4)
    n = 2
    vfacs = np.sqrt(np.array([1.5, 2, 2.5]))
    plt.figure()
    for vf in vfacs:
        data = pd.read_csv("data/v%.2fescvels2_%.3f_%.3f.txt"%(vf, 2, np.log10(nt)))
        r,v = mdimarr(data, nt, n)
        plt.plot(r[0,:,1], r[1,:,1], label = r"$v_{mult} = %.3f$"%vf)
    plt.legend()
    plt.axis("equal")
    plt.xlabel(r"$x$ [AU]")
    plt.ylabel(r"$y$ [AU]")
    plt.show()
    plt.figure()


def jupiterplots():
    nt = 20000
    n = 3
    massms = np.array([1,10, 1000])

    for i, title in enumerate(["Fixed sun", r"Non-fixed sun. $M_J * 1000$"]):
        plt.figure(i)
        plt.gca().set_aspect("equal", adjustable = 'box')
        plt.title(title)
        plt.xlabel("x[AU]")
        plt.ylabel("y[AU]")

    plt.figure(0)
    for i in range(3):
        fdata = pd.read_csv("data/m%.1fjupiter3_2.000_%.3f.txt"%(massms[i], np.log10(nt)))
        r,v = mdimarr(fdata, nt, n)
        plt.plot(r[0, :, 1], r[1, :, 1],
                 label = "%s, multiplier %i"%("earth",massms[i]))
    plt.plot(r[0, :, 2], r[1, :, 2], label = "Jupiter")
    plt.legend()
    plt.axis([-10, 10, -10, 10])
    #plt.savefig("figures/fsunjup.pdf")

    plt.figure(1)
    nfdata = pd.read_csv("data/nfm%.1fjupiter3_2.000_%.3f.txt"%(1000, np.log10(nt)))
    r,v = mdimarr(nfdata, nt, n)
    for i,pname in enumerate(['sun', 'earth', 'jupiter']):
        plt.plot(r[0, :, i], r[1, :, i], label = pname)
    plt.legend()
    plt.axis([-5, 5, -5, 5])
    #plt.savefig("figures/nfsunjup.pdf")

    energy_data = pd.read_csv("data/energy_nfjupiter3_2.000_%.3f.txt"%np.log10(nt))
    t = energy_data["      timestep"]
    PE = energy_data["PE"]
    KE = energy_data["KE"]
    plt.figure(2)
    plt.title("Total energies with high mass Jupiter, non-fixed sun")
    plt.plot(t, PE, label = r"$E_P$")
    plt.plot(t, KE, label = r"$E_K$")
    plt.plot(t, PE-KE, label = r"$E_P - E_K$")
    plt.xlabel("timestep")
    plt.yscale("log")
    plt.legend()
    plt.show()


def solarplots():
    nt = 80000
    n = 10
    data = pd.read_csv("data/solar%i_2.000_%.3f.txt"%(n, np.log10(nt)))
    r,v = mdimarr(data, nt, n)
    plt.figure()
    plt.title("Solar System on xy-plane")
    for i in range(n):
        plt.plot(r[0, :, i], r[1, :, i], label = "planet %i"%i)
    #plt.legend()
    plt.gca().set_aspect("equal", adjustable = 'box')
    plt.show()

    energy_data = pd.read_csv("data/energy_solar%i_2.000_%.3f.txt"%(n, np.log10(nt)))
    t = energy_data["      timestep"]
    PE = energy_data["PE"]
    KE = energy_data["KE"]
    plt.figure()
    plt.title("Total Energy for all planets")
    plt.plot(t, PE, label = r"$E_P$")
    plt.plot(t, KE, label = r"$E_K$")
    plt.plot(t, PE-KE, label = r"$E_P - E_K$")
    plt.xlabel("timestep")
    plt.yscale("log")
    plt.legend()
    plt.show()


def mercury_precession():
    """
    Plots mercury orbit and calculates precession?
    """
    nt = int(1e5)
    data = pd.read_csv("data/mercury2_2.000_%.3f.txt"%(np.log10(nt)))
    n = 2
    r, v = mdimarr(data, nt, n)
    r1 = r[:,:nt//100,1]
    r2 = r[:,99*nt//100:,1]
    plt.figure()
    plt.plot(0,0,".", color='orange', label='The Sun')
    plt.plot(r1[0], r1[1], ls ='-', color='k', label="1st year")
    plt.plot(r2[0], r2[1], ls ='dotted', color='lime', label="100th year")
    r_min1 = np.argmin(r1[0]**2 + r1[1]**2)
    x_p1 = r1[0, r_min1]
    y_p1 = r1[1, r_min1]
    r_min2 = np.argmin(r2[0]**2 + r2[1]**2)
    x_p2 = r2[0, r_min2]
    y_p2 = r2[1, r_min2]
    plt.plot(x_p1, y_p1, 'x', color='k', label="1st perihelion")
    plt.plot(x_p2, y_p2, '.', color='lime', label="100th perihelion")
    plt.xlabel(r"$x$ [AU]")
    plt.ylabel(r"$y$ [AU]")
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.axis("equal")
    plt.savefig("figures/mercury.pdf")
    plt.show()
    theta1 = np.arctan(y_p1/x_p1)
    theta2 = np.arctan(y_p2/x_p2)

    theta_p = (theta2 - theta1)*206265

    print("Precession angle: %.2f\"" %theta_p)

def plot_3D():
    n = 10
    nt = 80000
    data = pd.read_csv("data/solar10_2.000_4.903.txt")
    r,v = mdimarr(data, nt, n)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for j in range(n):
        ax.plot(r[0,:,j],r[1,:,j],r[2,:,j])

    ax.set_xlabel(r"$x$ [AU]")
    ax.set_ylabel(r"$y$ [AU]")
    ax.set_zlabel(r"$z$ [AU]")
    ax.set_zlim(-20,20)
    ax.set_zticks([-10,0,10])
    fig.savefig("figures/solar_system_3D.pdf")

#betaplots()
#mercury_precession()
#plot_3D()
#err_plots()

if __name__ == "__main__":
    #call functions here as wanted
    pass
