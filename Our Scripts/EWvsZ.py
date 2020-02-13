
import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math
from Redshift import redshift as z

z = np.linspace(0.051, 14, 1000)
#print (z)
EW=[mod.EW(i) for i in z]
#print (t_rec)

plotter.plot(z,EW,"A plot of redshift against EW","z","EW")