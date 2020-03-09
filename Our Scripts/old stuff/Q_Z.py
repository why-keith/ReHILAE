"""
Generates a plot of the evolution of ionised hydrgoen against redshift
"""
import numpy as np
from scipy.integrate import odeint
import MPhys_model as mod
import matplotlib.pyplot as plt
from Redshift import redshift as z
import plot_saver as save
import plotter
import random_array_generator as rag
import pandas as pd

if __name__=="__main__":
    path=None
else:
    path=save.folder

#TIME CONDITIONS
startT = mod.startT
finishT = mod.finishT
TStep = mod.TStep

#GENERATE Z AND t ARRAYS
z,t = z(startT, finishT, TStep)

#RETURNS dQ_dt 
def Q_Hii_dot(Q,t):
    dQ_dt = mod.Q_Hii_dot_UV(mod.z(t),Q)
    return dQ_dt

#GENERATE Q ARRAY
Q = odeint(Q_Hii_dot, mod.Q_Hii_zero, t)
Q[Q>1.0] = 1.0 # 100% HII
Q[Q<0.0] = 0.0 # 100% HI

#PLOT Z AGAINST Q
plt.plot(z,Q)
plt.show()