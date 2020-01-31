import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math
import Redshift as z

#TIME VALUES
startTime=0.01
finishTime=14
timeStep = 0.01

#GENERATE Z ARRAY
Z = z.redshift(startTime, finishTime, timeStep)

#GENERATE f_esc ARRAY
f = []
for i in Z:
    f.append(mod.f_esc(i))

#QUICK PLOT
plotter.plot(Z,f,"A plot of f_esc against redshift ","Redshift","f_esc")