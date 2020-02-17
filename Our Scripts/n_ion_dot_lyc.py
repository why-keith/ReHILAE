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
tStep = mod.tStep

#GENERATE Z AND t ARRAYS
Z,t = z(startT, finishT, tStep)

#GENERATE f_esc ARRAY
n_ion = []
for i in Z:
    n_ion.append(mod.n_ion_dot_LyC(i))

#PLOTS Z AGAINST f
plotter.plot(Z,n_ion,"A plot of the production rate of lyman continuum photons against redshift ","Redshift","n_ion_dot_lyc",path)
