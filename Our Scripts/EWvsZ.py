
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
z,t = z(startZ, finishZ, zStep)

EW=[mod.EW(i) for i in z]

plotter.plot(z,EW,"A plot of the equivalent width against redshift","Redshift","EW")