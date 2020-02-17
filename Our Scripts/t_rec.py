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

if __name__=="__main__":
    path=None
else:
    path=save.folder

#Z VALUES
startZ = mod.startZ
finishZ = mod.finishZ
zStep = mod.zStep

#GENERATE Z AND t ARRAYS
Z,t = z(startZ, finishZ, zStep)

t_rec=[math.log10(mod.t_rec(i)) for i in Z]
#print (t_rec)

plotter.plot(Z,t_rec,"A plot of recombination time against redshift","Redshift","log(t_rec) / log(s)",path)