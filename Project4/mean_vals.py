import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

from numpy import cumsum
from cycler import cycler

plt.rcParams.update({'font.size': 14})
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')

L = float(input("Dimension of spin matrix: "))
temp = float(input("Temperature: "))
tol = float(input("Randomness [0,1]: "))

data = pd.read_csv("data/mcdep_%d_%.2f_%.2f.csv"%(L,temp,tol))

E = data["E"]
mcs = data["mcs"]+1
M = data["M"]
acpt = data["acpt"]

N = len(acpt)-1

mE = cumsum(E)/mcs
mM = cumsum(M)/mcs

absM = cumsum(np.abs(M))/mcs

mE2 = cumsum(E**2)/mcs
mM2 = cumsum(M**2)/mcs

sigmaE2 = mE2 - mE**2

Cv = sigmaE2
chi = mM2 - absM**2

plt.figure("me")
plt.semilogx(mcs, mE/int(L**2))
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"$\left \langle E \right \rangle$")
plt.title("Mean energy per spin")
plt.tight_layout()
plt.savefig("figs/mean_energy_%d_%.2f_%.2f.pdf"%(L,temp,tol))
plt.close()

plt.figure("mm")
plt.semilogx(mcs, mM/int(L**2))
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"$\left \langle M \right \rangle$")
plt.title("Mean magnetisation per spin")
plt.tight_layout()
plt.savefig("figs/mean_magnetisation_%d_%.2f_%.2f.pdf"%(L,temp,tol))
plt.close()

plt.figure("m|m|")
plt.semilogx(mcs, absM/int(L**2))
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"$\left \langle |M| \right \rangle$")
plt.title("Absolute magnetisation per spin")
plt.tight_layout()
plt.savefig("figs/abs_magnetisation_%d_%.2f_%.2f.pdf"%(L,temp,tol))
plt.close()

plt.figure("me2")
plt.semilogx(mcs, mE2/int(L**2))
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"$\left \langle E^2 \right \rangle$")
plt.title("Second moment of energy")
plt.tight_layout()
plt.savefig("figs/second_E_%d_%.2f_%.2f.pdf"%(L,temp,tol))
plt.close()

plt.figure("mm2")
plt.semilogx(mcs, mM2/int(L**2))
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"$\left \langle M^2 \right \rangle$")
plt.title("Second moment of magnetisation")
plt.tight_layout()
plt.savefig("figs/second_M_%d_%.2f_%.2f.pdf"%(L,temp,tol))
plt.close()

plt.figure("cv")
plt.semilogx(mcs, Cv/int(L**2))
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"$\left \langle C_v \right \rangle$")
plt.title("Mean specific heat")
plt.tight_layout()
plt.savefig("figs/Cv_%d_%.2f_%.2f.pdf"%(L,temp,tol))
plt.close()

plt.figure("acpt")
plt.plot(mcs[:1000], acpt[:1000])
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"accepted")
plt.title("Configuration changes")
plt.tight_layout()
plt.savefig("figs/acpt_%d_%.2f_%.2f.pdf"%(L,temp,tol))
plt.close()


print("Monte Carlo cycles: %d %d %d" %(mcs[int(N/100)], mcs[int(N/10)], mcs[N]))
print("Mean energy per spin: %.6f %.6f %.6f" %(mE[int(N/100)]/int(L**2),
                                                      mE[int(N/10)]/int(L**2), mE[N]/int(L**2)))
print("Mean magnetisation per spin: %.6f %.6f %.6f" %(mM[int(N/100)]/int(L**2),
                                                      mM[int(N/10)]/int(L**2), mM[N]/int(L**2)))
print("Mean absolute magnetisation per spin: %.6f %.6f %.6f"
                            %(absM[int(N/100)]/int(L**2), absM[int(N/10)]/int(L**2), absM[N]/int(L**2)))
print("Mean specific heat per spin: %.7f %.7f %.7f" %(Cv[int(N/100)]/int(L**2),
                                                      Cv[int(N/10)]/int(L**2), Cv[N]/int(L**2)))
print("Magnetic susceptibility per spin: %.6e %.6e %.6e" %(chi[int(N/100)]/int(L**2),
                                                    chi[int(N/10)]/int(L**2), chi[N]/int(L**2)))
