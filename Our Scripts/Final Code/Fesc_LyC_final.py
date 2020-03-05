from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math
# DATA

def Function(x,a1,a2):
  """
    Function for the equation of the line that will be fitted to the data.
  """
  return a1*x + a2

xs = pylab.array([79, 129, 83, 98, 75, 29, 15, 4])                        # Equivalent Width
ys = pylab.array([0.132, 0.074, 0.072, 0.058, 0.056, 0.045, 0.032, 0.01]) # Escape Fraction
ys_err = 0.2*ys                                                           #Error in Escape Fraction

p2 = pylab.polyfit(xs, ys, 1.0)
p = pylab.poly1d(p2)
params = p.coefficients                                                   # initial estimation for parameters from polyfit
pfit, pcov = curve_fit(Function, xs,ys,p0=params,sigma=ys_err)
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
xdata = [i for i in range(0,160)]
data = [Function(i,a1,a2) for i in xdata]
plt.scatter(xs, ys, color='black')
plt.errorbar(xs, ys, yerr=ys_err, ls = 'none', color='black')
plt.plot(xdata, data, color='steelblue', label=r'$f_{esc,LyC}$ = (0.00064±0.00013)$EW_{0}$  + (0.00941±0.00364)')
plt.xlabel(r'EW [Å]')
plt.ylabel(r'$f_{esc,LyC}$')
plt.legend()
plt.show()