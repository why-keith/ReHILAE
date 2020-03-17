from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math

def Function(x, a1, a2):
    return a1*x + a2

x = pylab.array([2.5,2.8,2.9,3.1, 3.3, 3.7, 4.1, 4.5, 4.8, 5.0, 5.3, 5.8])
y = np.array([117.134349876, 122.113769531, 172.679885864, 143.5936203, 149.371124268, 96.3648490905, 189.002449036, 181.127731323, 95.0448436864, 125.935153962, 143.343048096, 251.93561])
y_err_up =pylab.array([186.35171, 160.96419, 321.31418, 382.54525, 258.18398, 70.23115, 763.57573, 223.72103, 100.98794, 195.63354, 159.8391, 620.82504])
y_err_down =pylab.array([53.82247, 53.55149, 90.78461, 72.25881, 74.65443, 39.08451, 96.45726, 85.27828, 42.09321, 52.08989, 71.00167, 186.34499])
y_err = [(y_err_up[i]-y_err_down[i])/2 for i in range(len(y_err_down))]
p2 = np.polyfit(x, y, 1.0)
p = pylab.poly1d(p2)
Params = p.coefficients
pfit, pcov = curve_fit(Function, x, y, p0=Params, sigma=y_err)
a1,a2= pfit[0], pfit[1]

error = []
for i in range(len(pfit)):
  try:
    error.append(np.absolute(pcov[i][i])**0.5)
  except:
    error.append(0.00)
BFP = pfit
err_BFP = np.array(error)

def EquivalentWidth(x):
    return 6.67298*x + 80.96867

def Fesc(z,a1,a2):
    a1*EquivalentWidth(z) + a2

print ('Function: y = (%s+-%s)x  + (%s+-%s)' %(round(BFP[0],5),round(err_BFP[0],5),round(BFP[1],5),round(err_BFP[1],5)))

a1, a2, = BFP[0], BFP[1]
a1_err, a2_err = err_BFP[0], err_BFP[1]
zs = list(np.linspace(0.0051,14,10000))
scipy_fit = [Fesc(z,a1,a2) for z in zs]

print(scipy_fit)

plt.plot(zs, scipy_fit)
plt.show()


