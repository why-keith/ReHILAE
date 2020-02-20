"""
Produces a plot of the production rate of lyman continuum photons against redshift
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

n_ion=[math.log10(mod.n_ion_dot_UV(i)) for i in z]

#PLOTS Z AGAINST f
plotter.plot(z,n_ion,"A plot of the production rate of lyman continuum photons against redshift ","Redshift (z)",r"$log(\dot{n}_{ion, LyC})[s^{-1}Mpc^{-3}]$",path)
