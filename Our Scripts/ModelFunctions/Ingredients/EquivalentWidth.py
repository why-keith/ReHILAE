from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math
import matplotlib

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


print ('Function: y = (%s+-%s)x  + (%s+-%s)' %(round(BFP[0],5),round(err_BFP[0],5),round(BFP[1],5),round(err_BFP[1],5)))

a1, a2, = BFP[0], BFP[1]
a1_err, a2_err = err_BFP[0], err_BFP[1]
zs = list(np.linspace(0.0051,20,10000))
scipy_fit = [Function(z,a1,a2) for z in zs]

n=1000
a1_list = np.random.normal(a1, a1_err*0.2, n)
a2_list = np.random.normal(a2, a2_err*0.2, n)

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

plt.figure('EW_Z')
#plt.plot(zs,scipy_fit,color='black',label='SciPy fit')
plt.errorbar(x,y,yerr=y_err,color='black',ls='none')
plt.fill_between(zs,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "grey", edgecolor = "black", linewidth = 1.2, label=r'$1\sigma$ conf. interval')
plt.plot(zs, median, label='ReHiLAE (this study, median)',color='black')
plt.scatter(x,y,color='black',marker='.',label='Sobral+2018')
plt.xlabel('Redshift (z)')
plt.ylabel(r'Ly$\alpha$ EW$_{0} \ [\mathring{A}$]')
plt.xlim(0,20)
plt.ylim(0)
plt.tick_params(which='both',direction='in',right=True,top=True)
plt.legend()
matplotlib.rcParams['lines.linewidth'] = 6
matplotlib.rcParams['axes.linewidth'] = 2.0
matplotlib.rcParams['xtick.major.size'] = 9
matplotlib.rcParams['xtick.minor.size'] = 5
matplotlib.rcParams['xtick.major.width'] = 1.9
matplotlib.rcParams['xtick.minor.width'] = 1.3
matplotlib.rcParams['ytick.major.size'] = 9
matplotlib.rcParams['ytick.minor.size'] = 4
matplotlib.rcParams['ytick.major.width'] = 1.9
matplotlib.rcParams['ytick.minor.width'] = 1.3
plt.show()