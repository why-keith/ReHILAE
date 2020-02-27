"""
Produces a plot of the availability of lyman continuum photons against redshift 
"""
import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math
from Redshift import redshift as z
import plot_saver as save
import random_array_generator as rag
import pandas as pd

if __name__=="__main__":
    path=None
else:
    path=save.folder

#TIME CONDITIONS
startT = mod.startT
finishT = mod.finishT
TStep = mod.TStep

#GENERATE Z AND t ARRAYS
z,t = z(startT, finishT, TStep)

#GENERATE f_esc ARRAY
Q = [math.log10(mod.Q_ion_LyC(i)) for i in z]

#PLOTS Z AGAINST f
plotter.plot(z,Q,"A plot of the availability of lyman continuum photons against redshift ","Redshift (z)",r"$Q_{ion, LyC}\times10^{50}$",path) 