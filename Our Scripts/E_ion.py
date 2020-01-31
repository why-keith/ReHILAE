import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math


z=np.arange(0,14,0.1)
#print (z)
E_ion=[math.log10(mod.E_ion(i)) for i in z]
#print (t_rec)

plotter.plot(z,E_ion,"A plot of redshift against ionisation efficiency","z","log(ξ_ion)")