from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math

def rhoLya(z):
    x = math.log10(1+z)
    if z < 5.8 :
      return 1.13869*x + 39.18853
    elif z > 5.8 :
      return  -5.288*x + 44.5388067

xdata = [i for i in range(0, 14)]
densities = [rhoLya(z) for z in xdata]

plt.plot(xdata, densities, color = 'steelblue')
plt.show()