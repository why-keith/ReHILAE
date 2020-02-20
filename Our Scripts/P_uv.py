import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math
from Redshift import redshift as z

#Z VALUES
startT = mod.startT
finishT = mod.finishT
TStep = mod.TStep

#GENERATE Z AND t ARRAYS
z,t = z(startT, finishT, TStep)

#GENERATE P_uv ARRAY
P_uv=[math.log10(mod.P_uv(i)) for i in z]

#PLOT Z AGAINST P_uv
plotter.plot(z,P_uv,"A plot of redshift against UV density","z",r"$\log(\rho_{UV})$")