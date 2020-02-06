import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math
from Redshift import redshift as z

#TIME VALUES
startTime = mod.startTime
finishTime = mod.finishTime
timeStep = mod.timeStep

#GENERATE Z AND t ARRAYS
Z,t = z(startTime, finishTime, timeStep)

#GENERATE P_uv ARRAY
P_uv=[mod.P_uv(i) for i in Z]

#PLOT Z AGAINST P_uv
plotter.plot(Z,P_uv,"A plot of redshift against UV density","z","log(P_uv)")