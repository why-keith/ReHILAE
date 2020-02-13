import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math
from Redshift import redshift as z

#Z VALUES
startZ = mod.startZ
finishZ = mod.finishZ
zStep = mod.zStep

#GENERATE Z AND t ARRAYS
Z,t = z(startZ, finishZ, zStep)


#GENERATE P_L_Lya ARRAY
P_L_Lya=[mod.P_L_Lya(i) for i in Z]


#PLOTS Z AGAINST P_L_Lya
plotter.plot(Z,P_L_Lya,"A plot of redshift against UV density","z","P_uv")