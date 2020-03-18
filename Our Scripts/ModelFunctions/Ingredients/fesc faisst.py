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
y=f0*(1+x)**alpha    
y_err_up=(f0+f0_err)*(1+x)**(alpha+alpha_err)
y_err_down=(f0-f0_err)*(1+x)**(alpha-alpha_err)

#x = pylab.array([2.5,2.8,2.9,3.1, 3.3, 3.7, 4.1, 4.5, 4.8, 5.0, 5.3, 5.8])
#y = np.array([117.134349876, 122.113769531, 172.679885864, 143.5936203, 149.371124268, 96.3648490905, 189.002449036, 181.127731323, 95.0448436864, 125.935153962, 143.343048096, 251.93561])
#y_err_up =pylab.array([186.35171, 160.96419, 321.31418, 382.54525, 258.18398, 70.23115, 763.57573, 223.72103, 100.98794, 195.63354, 159.8391, 620.82504])
#y_err_down =pylab.array([53.82247, 53.55149, 90.78461, 72.25881, 74.65443, 39.08451, 96.45726, 85.27828, 42.09321, 52.08989, 71.00167, 186.34499])
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


print ('Function: y = (%s+-%s)x  + (%s+-%s)' %(round(BFP[0],5),round(err_BFP[0],5),round(BFP[1],5),round(err_BFP[1],5)))

a1, a2, = BFP[0], BFP[1]
a1_err, a2_err = err_BFP[0], err_BFP[1]
zs = list(np.linspace(0.0051,14,10000))
scipy_fit = [Function(z,a1,a2) for z in zs]

a1_list = np.random.normal(a1, a1_err*0.2, 100)
a2_list = np.random.normal(a2, a2_err*0.2, 100)

all_runs = []
for _ in a1_list:
    all_runs.append([])

for A,B,loop in zip(a1_list,a2_list,all_runs):
    for z in zs:
        parameter = np.random.randint(2)
        if parameter == 0:
            loop.append(Function(z,A,a2))
        else:
            loop.append(Function(z,a1,B))


median, median_upper_percentile, median_lower_percentile = [],[],[]
for z in zs:
    ind = zs.index(z)
    equiWidths = [i[ind] for i in all_runs]
    median_lower_percentile.append(np.percentile(equiWidths,16))
    median.append(np.median(equiWidths))
    median_upper_percentile.append(np.percentile(equiWidths,84))


plt.figure(r'$f_{esc}$ - Faisst')
#plt.plot(zs,scipy_fit,color='black',label='SciPy fit')
#plt.scatter(x,y,color='black',marker='.')
#plt.errorbar(x,y,yerr=y_err,color='black',ls='none')
plt.fill_between(zs,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "grey", edgecolor = "black", linewidth = 1.2, label=r'68% Confidence Interval')
plt.plot(zs, median, "--", label='Median',color='blue')
plt.xlabel(r"$f_{esc}$")
plt.ylabel("Redshift")
plt.show()

