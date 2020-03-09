"""
Generates a plot of the escape fraction against redshift
"""
import matplotlib.pyplot as plt
import numpy as np
import Main_UV as mod
import math
import random_array_generator as rag
import pandas as pd

#Z VALUES
startT = mod.startT
finishT = mod.finishT
TStep = mod.TStep

#GENERATE Z AND t ARRAYS
z,t = mod.redshift(startT, finishT, TStep)

#GENERATE E_ion ARRAY
E_ion=[math.log10(mod.E_ion(i)) for i in z]

#PLOTS Z AGAINST xi
plt.figure("xi_ion")
plt.title("A plot of the UV production efficeincy against redshift")
plt.plot(z, E_ion, color = 'steelblue')
plt.xlabel("Redshift(z)")
plt.ylabel(r"$\log_{10}(\xi_{ion})$[Hz erg$^{-1}$]")
plt.show()