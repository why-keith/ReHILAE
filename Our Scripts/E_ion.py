import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math
from Redshift import redshift as z

startTime=0.01
finishTime=14
timeStep = 0.01

Z,t = z(startTime, finishTime, timeStep)
#print (z)
E_ion=[math.log10(mod.E_ion(i)) for i in Z]
#print (t_rec)

plotter.plot(Z,E_ion,"A plot of redshift against ionisation efficiency","z","log(ξ_ion)")