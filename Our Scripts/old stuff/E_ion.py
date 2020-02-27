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

#GENERATE E_ion ARRAY
E_ion=[math.log10(mod.E_ion(i)) for i in Z]

#PLOTS Z AGAINST E_ion
plotter.plot(Z,E_ion,"A plot of redshift against ionisation efficiency","z","log(ξ_ion)")