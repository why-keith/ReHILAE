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
n_ion = []
for i in Z:
    n_ion.append(mod.n_ion_dot_LyC(i))

#PLOTS Z AGAINST f
plotter.plot(Z,n_ion,"A plot of f_esc against redshift ","Redshift","n_ion_dot_lyc")
