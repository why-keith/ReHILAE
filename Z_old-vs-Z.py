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


#GENERATE Z ARRAY
#Z = z.redshift(startTime, finishTime, timeStep)
Z,t = z(startTime, finishTime, timeStep)

#GENERATE Z_alt ARRAY
#Z = z.redshift(startTime, finishTime, timeStep)
Z_alt,t = z(startTime, finishTime, timeStep, True)

#PLOTS Z AGAINST Z_alt
plotter.plot(Z,Z_alt,"A plot of old redshift against new redshift","Redshift","Old Redshift")