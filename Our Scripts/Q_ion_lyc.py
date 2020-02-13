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
Q = []
for i in Z:
    Q.append(mod.Q_ion_LyC(i))

#PLOTS Z AGAINST f
plotter.plot(Z,Q,"A plot of the availability of lyman continuum photons against redshift ","Redshift","Q_ion_LyC")