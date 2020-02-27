"""
Generates a plot of the equivalent width against redshift
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

EW=[mod.EW(i) for i in z]

plotter.plot(z,EW,"A plot of the equivalent width against redshift","Redshift","EW [Å]",path)