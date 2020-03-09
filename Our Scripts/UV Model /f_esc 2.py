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

#GENERATE f_esc ARRAY
f = [mod.f_esc_UV(i) for i in z]

#PLOTS Z AGAINST f
plt.figure("f_esc")
plt.title("A plot of the escape fraction against redshift")
plt.plot(z, f, color = 'steelblue')
plt.xlabel("Redshift(z)")
plt.ylabel(r"$f_{esc}$")
plt.show()