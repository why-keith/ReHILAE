"""
Produces a plot of the production rate of lyman continuum photons against redshift
"""
import matplotlib as plt
import numpy as np
import Main_LAE_Powerlaw as mod
import math
import matplotlib.pyplot as plt
import random_array_generator as rag
import pandas as pd


#TIME CONDITIONS
startT = mod.startT
finishT = mod.finishT
TStep = mod.TStep

#GENERATE Z AND t ARRAYS
z,t = mod.redshift(startT, finishT, TStep)

P1_1 = 1.13869
P1_2 = 39.18853
P2_1 = -5.288
P2_2 = 44.5388067
f1 = 0.00064
f2 = 0.00941

#GENERATE f_esc ARRAY
n_ion_1=[math.log(mod.n_ion_dot_LyC_PL_1(i, P1_1, P1_2, P2_1, P2_2, f1, f2)*2.938e+73) for i in z]
n_ion_2=[math.log(mod.n_ion_dot_LyC_PL_2(i, P1_1, P1_2, P2_1, P2_2, f1, f2)*2.938e+73) for i in z]


#PLOTS Z AGAINST f
plt.figure('n_ion_dot_LyC')
plt.title("A plot of the production rate of lyman continuum photons against redshift ")
plt.plot(z, n_ion_1, color = 'steelblue', label = 'powerlaw')
plt.plot(z, n_ion_2, color = 'black', label = 'powerlaw_corr')
plt.xlabel(r'Redshift (z)')
plt.ylabel(r"$log(\dot{n}_{ion, LyC})[s^{-1}Mpc^{-3}]$")
plt.legend()
plt.show()