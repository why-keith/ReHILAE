from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math

def Function(x, a1, a2):
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
a1_err, a2_err = err_BFP[0], err_BFP[1]

def EW(x):
    return 6.67298*x + 80.96867

def Fesc(z,a1,a2):
    return a1*EW(z) + a2

zs = list(np.linspace(0.0051,14,10000))
scipy_fit = [Fesc(z,a1,a2) for z in zs]

number_of_iterations = 100
a1_list = np.random.normal(a1, a1_err*0.2, number_of_iterations)
a2_list = np.random.normal(a2, a2_err*0.2, number_of_iterations)

all_runs = []
for _ in a1_list:
    all_runs.append([])

for A,B,loop in zip(a1_list,a2_list,all_runs):
    for z in zs:
        parameter = np.random.randint(2)
        if parameter == 0:
            loop.append(Fesc(z,A,a2))
        else:
            loop.append(Fesc(z,a1,B))


median, median_upper_percentile, median_lower_percentile = [],[],[]
for z in zs:
    ind = zs.index(z)
    equiWidths = [i[ind] for i in all_runs]
    median_lower_percentile.append(np.percentile(equiWidths,56))
    median.append(np.median(equiWidths))
    median_upper_percentile.append(np.percentile(equiWidths,34))

plt.figure('Fesc_LyC')
plt.plot(zs, scipy_fit, color='black', label='SciPy fit')
plt.fill_between(zs,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "grey", edgecolor = "black", linewidth = 1.2, label=r'68% Confidence Interval')
plt.plot(zs, median, "--", label='Median',color='blue')
plt.xlabel('Redshift (z)')
plt.ylabel(r'$f_{esc,LyC}$')
plt.legend()
plt.show()