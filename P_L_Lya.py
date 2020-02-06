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
print(Z)


#GENERATE P_L_Lya ARRAY
P_L_Lya=[mod.P_L_Lya(i) for i in Z]


#PLOTS Z AGAINST P_L_Lya
plotter.plot(Z,P_L_Lya,"A plot of redshift against UV density","z","log(P_uv)")