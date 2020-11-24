"""
Script for calculating thermodynamical quantities for 2D ising model
"""
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

from numpy import cumsum
from cycler import cycler

warnings.filterwarnings("ignore")

plt.rcParams.update({'font.size': 14})
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')

L = float(input("Dimension of spin matrix: "))
temp = float(input("Temperature: "))
tol = float(input("Configuration (0 or 0.5): "))

random_ordered = input("Add another configuration: [Y/n]: ")
test = False
data2 = False

data = pd.read_csv("data/mcdep_%d_%.2f_%.2f.csv"%(L, temp, tol))

E = data["E"]
mcs = data["mcs"]+1
M = data["M"]
acpt = data["acpt"]

N = len(mcs)-1

mE = cumsum(E)/mcs
mM = cumsum(M)/mcs

absM = cumsum(np.abs(M))/mcs

mE2 = cumsum(E**2)/mcs
mM2 = cumsum(M**2)/mcs

sigmaE2 = mE2 - mE**2
sigmaM2 = mM2 - absM**2

Cv = sigmaE2/temp**2
chi = sigmaM2/temp**2

if random_ordered in ('Y', 'y'):
    test = True
    tol2 = float(input("Configuration (0 or 0.5): "))
    data2 = pd.read_csv("data/mcdep_%d_%.2f_%.2f.csv"%(L, temp, tol2))

    E_2 = data2["E"]
    M_2 = data2["M"]
    acpt_2 = data2["acpt"]

    mE_2 = cumsum(E_2)/mcs
    mM_2 = cumsum(M_2)/mcs

    absM_2 = cumsum(np.abs(M_2))/mcs

    mE2_2 = cumsum(E_2**2)/mcs
    mM2_2 = cumsum(M_2**2)/mcs

    sigmaE2_2 = mE2_2 - mE_2**2

    Cv_2 = sigmaE2_2
    chi_2 = mM2_2 - absM_2**2


plt.figure("me")
fig, ax1 = plt.subplots()
ax1.semilogx(mcs, mE/int(L**2), label='random', color='g')
ax1.semilogx(mcs, mE/int(L**2), label='ordered', color='b')
ax1.set_xlabel("Monte carlo cycles")
ax1.set_ylabel(r"$\left \langle E \right \rangle$")
ax1.set_title("Mean energy per spin")
ax1.tick_params(axis='y', labelcolor='b')
if test:
    ax2 = ax1.twinx()
    ax2.semilogx(mcs, mE_2/int(L**2), color='g')
    ax2.tick_params(axis='y', labelcolor='g')
    ax1.legend(loc="center right")
plt.tight_layout()
plt.savefig("figs/mean_energy_%d_%.2f_%.2f.pdf"%(L, temp, tol))
plt.close()

plt.figure("mm")
plt.semilogx(mcs, mM/int(L**2), label='ordered')
if test:
    plt.semilogx(mcs, mM_2/int(L**2), label='random')
    plt.legend()
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"$\left \langle M \right \rangle$")
plt.title("Mean magnetisation per spin")
plt.tight_layout()
plt.savefig("figs/mean_magnetisation_%d_%.2f_%.2f.pdf"%(L, temp, tol))
plt.close()

plt.figure("m|m|")
plt.semilogx(mcs, absM/int(L**2), label='ordered')
if test:
    plt.semilogx(mcs, absM_2/int(L**2), label='random')
    plt.legend()
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"$\left \langle |M| \right \rangle$")
plt.title("Absolute magnetisation per spin")
plt.tight_layout()
plt.savefig("figs/abs_magnetisation_%d_%.2f_%.2f.pdf"%(L, temp, tol))
plt.close()

plt.figure("sigmaE")
plt.loglog(mcs, sigmaE2/int(L**2), label='ordered')
if test:
    plt.loglog(mcs, sigmaE2_2/int(L**2), label='random')
    plt.legend()
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"$\left \langle \sigma_E^2 \right \rangle$")
plt.title("Variance in energy")
plt.tight_layout()
plt.savefig("figs/sigmaE_%d_%.2f_%.2f.pdf"%(L, temp, tol))
plt.close()

plt.figure("cv")
plt.semilogx(mcs, Cv/int(L**2), label='ordered')
if test:
    plt.semilogx(mcs, Cv_2/int(L**2), label='random')
    plt.legend()
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"$\left \langle C_v \right \rangle$")
plt.title("Mean specific heat")
plt.tight_layout()
plt.savefig("figs/Cv_%d_%.2f_%.2f.pdf"%(L, temp, tol))
plt.close()

plt.figure("acpt")
plt.loglog(mcs, acpt, label='ordered')
if test:
    plt.loglog(mcs, acpt_2, label='random')
    plt.legend()
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"accepted")
plt.title("Configuration changes")
plt.tight_layout()
plt.savefig("figs/acpt_%d_%.2f_%.2f.pdf"%(L, temp, tol))
plt.close()

plt.figure("prob")
counts, bins = np.histogram(E[200:]/L**2, bins=int((E[200:].max()
                                                    -E[200:].min())/4)+1)
plt.hist(bins[:-1], bins, weights=counts/(N-199))
if temp > 2:
    plt.axvline(mE[N]/int(L**2)+np.sqrt(sigmaE2[N])/L**2, color="r", ls="--")
    plt.axvline(mE[N]/int(L**2)-np.sqrt(sigmaE2[N])/L**2, color="r", ls="--")
plt.xlabel("Energy per spin")
plt.ylabel("Probability")
plt.title("Probability distribution")
plt.tight_layout()
plt.savefig("figs/prob_%d_%.2f_%.2f.pdf"%(L, temp, tol))
plt.close()

print("Monte Carlo cycles: %d %d %d" %(mcs[int(N/100)], mcs[int(N/10)], mcs[N]))
print("Mean energy per spin: %.6f %.6f %.6f"
      %(mE[int(N/100)]/int(L**2), mE[int(N/10)]/int(L**2), mE[N]/int(L**2)))
print("Mean magnetisation per spin: %.6f %.6f %.6f"
      %(mM[int(N/100)]/int(L**2), mM[int(N/10)]/int(L**2), mM[N]/int(L**2)))
print("Mean absolute magnetisation per spin: %.6f %.6f %.6f"
      %(absM[int(N/100)]/int(L**2), absM[int(N/10)]/int(L**2),
        absM[N]/int(L**2)))
print("Mean specific heat per spin: %.7f %.7f %.7f"
      %(Cv[int(N/100)]/int(L**2), Cv[int(N/10)]/int(L**2), Cv[N]/int(L**2)))
print("Magnetic susceptibility per spin: %.6e %.6e %.6e"
      %(chi[int(N/100)]/int(L**2), chi[int(N/10)]/int(L**2), chi[N]/int(L**2)))
print("Standard deviation in energy per spin: %.6f %.6f %.6f"
      %(np.sqrt(sigmaE2[int(N/100)])/L**2, np.sqrt(sigmaE2[int(N/10)])/L**2,
        np.sqrt(sigmaE2[N])/L**2))
