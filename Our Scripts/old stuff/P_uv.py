import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math
from Redshift import redshift as z

#TIME VALUES
startZ = mod.startZ
finishZ = mod.finishZ
zStep = mod.zStep

#GENERATE Z AND t ARRAYS
Z,t = z(startZ, finishZ, zStep)

#GENERATE P_uv ARRAY
P_uv=[mod.P_uv(i) for i in Z]

#PLOT Z AGAINST P_uv
plotter.plot(Z,P_uv,"A plot of redshift against UV density","z","log(P_uv)")