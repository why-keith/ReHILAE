"""
Produces a plot of the production rate of lyman continuum photons against redshift
"""
import matplotlib as plt
import numpy as np
import Main_LAE as mod
import plotter
import math
import matplotlib.pyplot as plt
from Redshift import redshift as z
import plot_saver as save
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

P1_1 = -0.48111
P1_2 = 40.10956

f1 = 0.00064
f1_error = [0.00013]

f2 = 0.00941
f2_error = [0.00364]

P2_1 = -0.04764
P2_1_error = 0.00254

P2_2 = 0.45297
P2_2_error = [0.03671]

P2_3 = 38.97568
P2_3_error = [0.08654]

 

#GENERATE f_esc ARRAY
n_ion=[math.log(mod.n_ion_dot_LyC_PL(i, P1_1, P1_2, f1, f2)*2.938e+73) for i in z]

n_ion_Q = [math.log(mod.n_ion_dot_LyC_Q(i, P2_1, P2_2, P2_3, f1, f2)*2.938e+73) for i in z]

#PLOTS Z AGAINST f
plt.plot(z, n_ion, color = 'steelblue', label = 'powerlaw')
plt.plot(z, n_ion_Q, color = 'black', label = 'Quadratic')
plt.legend()
plt.show()