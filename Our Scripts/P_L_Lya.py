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
TStep = mod.TStep

#GENERATE Z AND t ARRAYS
z,t = z(startT, finishT, TStep)

#GENERATE P_L_Lya ARRAY
P_L_Lya=[math.log10(mod.P_L_Lya(i)) for i in z]


#PLOTS Z AGAINST P_L_Lya
plotter.plot(z,P_L_Lya,"A plot of the luminosity density of lyman alpha against redshift","Redshift (z)",r"$\rho_{L_{Lya}} [erg s^{−1} Mpc^{−3}]$",path)
