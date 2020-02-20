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

#GENERATE E_ion ARRAY
E_ion=[math.log10(mod.E_ion(i)) for i in z]

#PLOTS Z AGAINST E_ion
plotter.plot(z,E_ion,"A plot of redshift against ionisation efficiency","Redshift (z)",r"$\log(\xi_{ion})$[Hz erg$^{-1}$]")