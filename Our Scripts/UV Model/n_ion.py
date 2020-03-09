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

#GENERATE f_esc ARRAY
n_ion=[math.log10(mod.n_ion_dot_UV(i, P1, P2)*2.938e+73) for i in z]


#PLOTS Z AGAINST f
plt.figure("n_ion")
plt.title("A plot of the production rate of lyman continuum photons against redshift ")
plt.plot(z, n_ion, color = 'steelblue')
plt.xlabel("Redshift(z)")
plt.ylabel(r"$log(\dot{n}_{ion})[s^{-1}Mpc^{-3}]$")
plt.show()