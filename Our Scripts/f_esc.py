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

#GENERATE f_esc ARRAY
f = []
for i in Z:
    f.append(mod.f_esc(i))

#PLOTS Z AGAINST f
plotter.plot(Z,f,"A plot of f_esc against redshift ","Redshift","f_esc")
