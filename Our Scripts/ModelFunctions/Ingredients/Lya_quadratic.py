from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math

def Function(x,a1,a2,a3):
	return a1*x**2.+a2*x+a3

xs = [2.2,2.5,2.8,3.0,3.2, 3.3, 3.7, 4.1, 4.6, 4.8, 5.1, 5.3, 5.8 ]
y = np.array([0.52, 0.74, 0.77, 0.88, 0.84, 0.85, 1.01, 0.87, 1.19, 1.12, 1.27, 1.08, 1.10]) # data from SC4K Sobral 
ys = np.array([math.log10(i*10**40) for i in y])
ys_err = [0.05,0.075,0.095,0.095,0.09,0.095,0.18,0.13,0.35,0.32,0.4,0.185,0.185]
p2 = pylab.polyfit(xs, ys, 2.0)
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

print ('Power Law: y = (%s+-%s)x  + (%s+-%s)' %(round(BFP[0],5),round(err_BFP[0],5),round(BFP[1],5),round(err_BFP[1],5)))

a1, a2, a3 = BFP[0], BFP[1], BFP[2]
a1_err, a2_err,  a3_err = err_BFP[0], err_BFP[1], err_BFP[2]
zs = list(np.linspace(0.0051,14,10000))
scipy_fit = [Function(z,a1,a2,a3) for z in zs]

a1_list = np.random.normal(a1, a1_err*0.2, 500)
a2_list = np.random.normal(a2, a2_err*0.2, 500)
a3_list = np.random.normal(a3, a3_err*0.2, 500)

all_runs = []
for _ in a1_list:
    all_runs.append([])

for A,B,C,loop in zip(a1_list,a2_list,a3_list,all_runs):
    for z in zs:
        parameter = np.random.randint(3)
        if parameter == 0:
            loop.append(Function(z,A,a2,a3))
        elif parameter == 1:
            loop.append(Function(z,a1,B,a3))
        else:
            loop.append(Function(z,a1,a2,C))


median, median_upper_percentile, median_lower_percentile = [],[],[]
for z in zs:
    ind = zs.index(z)
    P_Lya = [i[ind] for i in all_runs]
    median_lower_percentile.append(np.percentile(P_Lya,16))
    median.append(np.median(P_Lya))
    median_upper_percentile.append(np.percentile(P_Lya,84))


plt.figure('Lya_Quad')
plt.plot(zs,scipy_fit,color='black',label='SciPy fit')
plt.scatter(xs,ys,color='black',marker='.')
plt.errorbar(xs,ys,yerr=ys_err,color='black',ls='none')
plt.fill_between(zs,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "grey", edgecolor = "black", linewidth = 1.2, label=r'68% Confidence Interval')
plt.plot(zs, median, "--", label='Median',color='blue')
plt.legend()
plt.show()