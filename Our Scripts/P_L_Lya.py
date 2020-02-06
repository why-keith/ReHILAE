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

P_L_Lya=[mod.P_L_Lya(i) for i in Z]

plotter.plot(Z,P_L_Lya,"A plot of redshift against Luminosity density of Lyman alpha","z","P_L_Lya")