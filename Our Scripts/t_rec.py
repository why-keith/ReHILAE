"""
Generates a plot of recombination time against redshift
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

t_rec=[math.log10(mod.t_rec(i)) for i in z]
#print (t_rec)

plotter.plot(z,t_rec,"A plot of recombination time against redshift","Redshift (z)",r"$\log_{10}(t_{rec}[s^{-1}])$",path)
