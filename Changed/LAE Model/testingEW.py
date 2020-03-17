from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math

def Function(x,a1,a2):
  """
    Function for the equation of the line that will be fitted to the data.
  """
  return a1*x + a2

xs = pylab.array([79, 129, 83, 98, 75, 29, 15, 4])                       
ys = pylab.array([0.132, 0.074, 0.072, 0.058, 0.056, 0.045, 0.032, 0.01]) 
ys_err = 0.2*ys                                                           
p2 = pylab.polyfit(xs, ys, 1.0)
p = pylab.poly1d(p2)
params = p.coefficients
pfit, pcov = curve_fit(Function, xs,ys,p0=params,sigma=ys_err)
error = [] 
for i in range(len(pfit)):
  try:
    error.append(np.absolute(pcov[i][i])**0.5)
  except:
    error.append(0.00)
BFP = pfit  
err_BFP = np.array(error) 

print ('Function: y = (%s+-%s)x  + (%s+-%s)' %(round(BFP[0],5),round(err_BFP[0],5),round(BFP[1],5),round(err_BFP[1],5)))

a1, a2, = BFP[0], BFP[1]
a1_error, a2_error = err_BFP[0], err_BFP[1]
xdata = list(np.linspace(0,200,10000))
data = [Function(i,a1,a2) for i in xdata]

a1_list = np.random.normal(a1, a1_error, 500)
a2_list = np.random.normal(a2, a2_error, 500)

all_runs = []
for _ in a1_list:
    all_runs.append([])

for A,B,loop in zip(a1_list,a2_list,all_runs):
    for EW in xdata:
        loop.append(Function(EW,A,B))

median, median_upper_percentile, median_lower_percentile = [],[],[]
for x in xdata:
    ind = xdata.index(x)
    f_escs = [i[ind] for i in all_runs]
    median_lower_percentile.append(np.percentile(f_escs,16))
    median.append(np.median(f_escs))
    median_upper_percentile.append(np.percentile(f_escs,84))

err_up = Function(148.9705, a1+err_BFP[0], a2+err_BFP[1])
err_down = Function(148.9705, a1-err_BFP[0], a2-err_BFP[1])
err = [(err_up-err_down)*0.5]

plt.figure("Fesc_LyC")
plt.scatter(xs, ys, color='black')
plt.errorbar(xs, ys, yerr=ys_err, ls = 'none', color='black')
plt.plot(xdata, data, color='black', label='SciPy fit')

plt.fill_between(xdata,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "grey", edgecolor = "black", linewidth = 1.2, label=r'$1\sigma$')
plt.scatter(148.9705, Function(148.9705,a1,a2), color='red')
plt.errorbar([148.9705], [Function(148.9705,a1,a2)], yerr=err , color='red', ls='none')
plt.plot(xdata, median, "--", label='median',color='blue')
plt.xlabel(r'Ly$\alpha$ EW [Å]')
plt.ylabel(r'$f_{esc,LyC}$')
plt.legend()
plt.show()