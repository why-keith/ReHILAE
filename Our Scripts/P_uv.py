import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math

z=np.arange(0,14,0.1)
#print (z)
P_uv=[mod.P_uv(i) for i in z]
#print (t_rec)

plotter.plot(z,P_uv,"A plot of redshift against UV density","z","log(P_uv)")