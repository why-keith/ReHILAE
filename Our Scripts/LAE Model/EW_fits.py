from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math


"""
SC4K data shown in Calhau et al. 2019 Table C3

available at: https://lancaster.app.box.com/s/t75t3v713yuibkvjk3ioqdesunxpv1fb/file/618125008170
"""

def linear(x, a1, a2):
    return a1*x + a2

x = pylab.array([2.5,2.8,2.9,3.1, 3.3, 3.7, 4.1, 4.5, 4.8, 5.0, 5.3, 5.8])
y = np.array([117.134349876, 122.113769531, 172.679885864, 143.5936203, 149.371124268, 96.3648490905, 189.002449036, 181.127731323, 95.0448436864, 125.935153962, 143.343048096, 251.93561])
y_err_up =pylab.array([186.35171, 160.96419, 321.31418, 382.54525, 258.18398, 70.23115, 763.57573, 223.72103, 100.98794, 195.63354, 159.8391, 620.82504])
y_err_down =pylab.array([53.82247, 53.55149, 90.78461, 72.25881, 74.65443, 39.08451, 96.45726, 85.27828, 42.09321, 52.08989, 71.00167, 186.34499])
y_err = [(y_err_up[i]-y_err_down[i])/2 for i in range(len(y_err_down))]
p2 = np.polyfit(x, y, 1.0)
p = pylab.poly1d(p2)
Params = p.coefficients
pfit, pcov = curve_fit(linear, x, y, p0=Params, sigma=y_err)
a1,a2= pfit[0], pfit[1]

error = [] # Empty error array
for i in range(len(pfit)): # Let's sample each individual parameter used
  try:
    error.append(np.absolute(pcov[i][i])**0.5)
  except:
    error.append(0.00)
BFP = pfit  # Best fit parameters
err_BFP = np.array(error) # Errors on best fit parameters

###### Print out the results for the best fit
print ('Function: y = (%s+-%s)x  + (%s+-%s)' %(round(BFP[0],5),round(err_BFP[0],5),round(BFP[1],5),round(err_BFP[1],5)))
# Show these on the plot as well

a1, a2, = BFP[0], BFP[1]
xdata = [i for i in range(0,16)]
data2 = [np.mean(y) for i in xdata]
data = [linear(i,a1,a2) for i in xdata]

plt.figure("EW_fits")
plt.scatter(x, y, color='black')
plt.errorbar(x, y, yerr=y_err, ls = 'none', color='black')
plt.plot(xdata, data2, color = 'mediumblue', label = ('Mean EW = %s' %(round(np.mean(y), 2))))
plt.plot(xdata, data, color='steelblue', label=r'EW Linear fit')
plt.xlabel(r'Redshift (z)')
plt.ylabel('EW [Å]')
plt.legend()
plt.show()