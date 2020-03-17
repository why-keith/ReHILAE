from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math

def Function(x, a1, a2):
    return a1*x + a2

f0=2.3/100
f0_err=0.1/100
alpha=1.17
alpha_err=0.02

x=np.linspace(0.0051,14,10000)
y=np.array([24.4+math.log10(1+i) for i in x])   

plt.figure()
#plt.plot(zs,scipy_fit,color='black',label='SciPy fit')
#plt.scatter(x,y,color='black',marker='.')
#plt.errorbar(x,y,yerr=y_err,color='black',ls='none')
#plt.fill_between(zs,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "grey", edgecolor = "black", linewidth = 1.2, label=r'$1\sigma$')
plt.plot(x, y, "-", label='Median',color='blue')
plt.xlabel(r"$log(\xi_{ion}) / log(Hz erg^{-1}$")
plt.ylabel("Redshift")
plt.show()


