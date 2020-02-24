"""
Generates a plot of the escape fraction against redshift
"""
import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math
from Redshift import redshift as z
import plot_saver as save

if __name__=="__main__":
    path=None
else:
    path=save.folder

#Z VALUES
startT = mod.startT
finishT = mod.finishT
TStep = mod.TStep

#GENERATE Z AND t ARRAYS
z,t = z(startT, finishT, TStep)

#GENERATE f_esc ARRAY
f = []
for i in z:
    f.append(mod.f_esc_LyC(i))

#PLOTS Z AGAINST f
plotter.plot(z,f,"A plot of the escape fraction against redshift","Redshift (z)",r"$f_{esc}$",path)
