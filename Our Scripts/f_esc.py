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
EW = [i for i in range(250)]
f = [mod.f_esc_LyC(i) for i in EW]

#PLOTS Z AGAINST f
plotter.plot(EW,f,"A plot of the escape fraction against equivalent width","EW",r"$f_{esc}$",path)
