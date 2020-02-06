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

#GENERATE f_esc ARRAY
f = []
for i in Z:
    f.append(mod.f_esc(i))

#PLOTS Z AGAINST f
plotter.plot(Z,f,"A plot of f_esc against redshift ","Redshift","f_esc")
