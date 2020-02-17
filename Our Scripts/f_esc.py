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
tStep = mod.tStep

#GENERATE Z AND t ARRAYS
Z,t = z(startT, finishT, tStep)

#GENERATE f_esc ARRAY
f = []
for i in Z:
    f.append(mod.f_esc_UV(i))

#PLOTS Z AGAINST f
plotter.plot(Z,f,"A plot of the escape fraction against redshift","Redshift","f_esc",path)
