"""
Generates a plot of the luminosity density of lyman alpha against redshift
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


#GENERATE P_L_Lya ARRAY
P_L_Lya=[mod.P_L_Lya(i) for i in Z]


#PLOTS Z AGAINST P_L_Lya
plotter.plot(Z,P_L_Lya,"A plot of the luminosity density of lyman alpha against redshift","z","P_L_Lya",path)
