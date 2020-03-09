"""
Generates a plot of the escape fraction against redshift
"""
import matplotlib.pyplot as plt
import numpy as np
import Main_UV as mod
import math
import random_array_generator as rag
import pandas as pd

#TIME CONDITIONS
startT = mod.startT
finishT = mod.finishT
TStep = mod.TStep

#GENERATE Z AND t ARRAYS
z,t = mod.redshift(startT, finishT, TStep)
P1 = -5.288
P2 = 30.39894

#GENERATE P_uv ARRA
P=[mod.P_UV(i, P1, P2) for i in z]



#PLOTS Z AGAINST f
plt.figure("P_UV")
plt.title("A plot of redshift against UV density")
plt.plot(z, P, color = 'steelblue')
plt.xlabel("Redshift(z)")
plt.ylabel(r"$\log_{10}(\rho_{UV}) [erg s^{−1} Mpc^{−3}]$")
plt.show()