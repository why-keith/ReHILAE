import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math
from Redshift import redshift as z


##########################################
#Broke for now but like do we really care 
##########################################

#Z VALUES
startZ = mod.startZ
finishZ = mod.finishZ
zStep = mod.zStep

#GENERATE Z ARRAY
#Z = z.redshift(startTime, finishTime, timeStep)
Z,t = z(startZ, finishZ, zStep)

#GENERATE Z_alt ARRAY
#Z = z.redshift(startTime, finishTime, timeStep)
Z_alt,t = Z,t = z(startZ, finishZ, zStep, True)

#PLOTS Z AGAINST Z_alt
plotter.plot(Z,Z_alt,"A plot of old redshift against new redshift","Redshift","Old Redshift")