import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math
from Redshift import redshift as z

#TIME VALUES
startTime=mod.startTime
finishTime=mod.finishTime
timeStep = mod.timeStep

#GENERATE Z AND t ARRAYS
Z,t = z(startTime, finishTime, timeStep)

t_rec=[math.log10(mod.t_rec(i)) for i in Z]
#print (t_rec)

plotter.plot(Z,t_rec,"A plot of redshift against recombination time","z","log(t_rec) / log(s)")