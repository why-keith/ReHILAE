from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math
import random_array_generator as rag
# DATA


################################################################################
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
################################################################################
#iteration data for fesc
f1 = [0.00064]
f1_error = [0.00013]

f2 = [0.00941]
f2_error = [0.00364]
 
f1_f_esc = rag.random_Arrays(len(f1),f1,f1_error,f1_error)
f2_f_esc = rag.random_Arrays(len(f2),f2,f2_error,f2_error)

f_esc_iterations = []
for i in range(len(f1_f_esc)):
    dummy = []
    for j in xdata:
        dummy.append(f1_f_esc[i][0]*j + f2_f_esc[i][0])
    f_esc_iterations.append(dummy)
        
median, median_lower_percentile, median_upper_percentile = rag.median_y_values(len(f_esc_iterations[0]),f_esc_iterations)
################################################################################
plt.figure("Fesc_LyC")
plt.scatter(xs, ys, color='black')
plt.errorbar(xs, ys, yerr=ys_err, ls = 'none', color='black')
plt.plot(xdata, data, color='steelblue', label='')

plt.fill_between(xdata,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "grey", edgecolor = "black", linewidth = 1.2)
plt.scatter(148.9705, Function(148.9705,a1,a2), color='red')
err_up = Function(148.9705, a1+err_BFP[0], a2+err_BFP[1])
err_down = Function(148.9705, a1-err_BFP[0], a2-err_BFP[1])
err = [(err_up-err_down)*0.5]
plt.errorbar([148.9705], [Function(148.9705,a1,a2)], yerr=err , color='red', ls='none')
plt.plot(xdata, median, "--", label='median',color='black')
plt.xlabel(r'Ly$\alpha$ EW [Å]')
plt.ylabel(r'$f_{esc,LyC}$')
plt.legend()
plt.show()