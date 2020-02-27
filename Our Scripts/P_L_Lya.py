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

#GENERATE P_L_Lya ARRAY
P_L_Lya=[mod.P_L_Lya(i) for i in z]

#PLOTS Z AGAINST P_L_Lya
plotter.plot(z,P_L_Lya,"A plot of the luminosity density of lyman alpha against redshift","z","P_L_Lya",path)
