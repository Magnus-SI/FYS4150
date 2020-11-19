import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.integrate import cumtrapz

plt.rcParams.update({'font.size': 14})

data = pd.read_csv("data/mcdep_2_2.40_0.50.csv")

mcs = data["mcs"]
E = data["E"]
M = data["M"]
acpt = data["acpt"]

mE = cumtrapz(y=E, x=mcs, initial=E[0])/mcs
mM = cumtrapz(y=M, x=mcs, initial=M[0])/mcs

mE2 = cumtrapz(y=E**2, x=mcs, initial=(E**2)[0])/mcs
mM2 = cumtrapz(y=M**2, x=mcs, initial=(M**2)[0])/mcs

plt.semilogx(mcs, mE)
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"$\left \langle E \right \rangle$")
plt.title("Mean energy per spin")
plt.tight_layout()
plt.savefig("figs/mean_energy.pdf")
plt.show()

plt.semilogx(mcs, mM)
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"$\left \langle M \right \rangle$")
plt.title("Mean magnetisation per spin")
plt.tight_layout()
plt.savefig("figs/mean_magnetisation.pdf")
plt.show()

plt.semilogx(mcs, mE2)
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"$\left \langle E^2 \right \rangle$")
plt.title("Second moment of energy")
plt.tight_layout()
plt.savefig("figs/second_E.pdf")
plt.show()

plt.semilogx(mcs, mM2)
plt.xlabel("Monte carlo cycles")
plt.ylabel(r"$\left \langle M^2 \right \rangle$")
plt.title("Second moment of magnetisation")
plt.tight_layout()
plt.savefig("figs/second_M.pdf")
plt.show()
