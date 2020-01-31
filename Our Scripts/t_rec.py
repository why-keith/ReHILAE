import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math

z=np.arange(0,14,0.1)
#print (z)
t_rec=[math.log10(mod.t_rec(i)) for i in z]
#print (t_rec)

plotter.plot(z,t_rec,"A plot of redshift against recombination time","z","log(t_rec) / log(s)")