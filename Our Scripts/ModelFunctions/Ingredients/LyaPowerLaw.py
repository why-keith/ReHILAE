from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math
import matplotlib

def Function(x,a1,a2):
	return a1*x + a2

def rhoLya(z,a1,a2):
    x = math.log10(1+z)
    return a1*x + a2

xs = pylab.array([2.2,2.5,2.8,3.0,3.2, 3.3, 3.7, 4.1, 4.6, 4.8, 5.1, 5.3, 5.8])                         
y1 = pylab.array([0.52, 0.74, 0.77, 0.88, 0.84, 0.85, 1.01, 0.87, 1.19, 1.12, 1.27, 1.08, 1.10])
ys = [math.log10(i*10**40) for i in y1]
xs_log = pylab.array([math.log10(i+1) for i in xs])
y_sigma_up = [0.05, 0.08, 0.1, 0.1, 0.9, 0.1, 0.2, 0.15, 0.4, 0.37, 0.48, 0.21, 0.21]
y_sigma_down =[0.05, 0.07, 0.09, 0.09, 0.09, 0.09, 0.16, 0.11, 0.3, 0.27, 0.32, 0.16, 0.16]
y_sigma = [(y_sigma_up[i]+y_sigma_down[i])/2 for i in range(len(y_sigma_down))]
p2 = pylab.polyfit(xs_log, ys, 1.0)
p = pylab.poly1d(p2)
params = p.coefficients

pfit, pcov = curve_fit(Function, xs_log,ys,p0=params,sigma=y_sigma)

error = []
for i in range(len(pfit)):
  try:
    error.append(np.absolute(pcov[i][i])**0.5)
  except:
    error.append(0.00)
BFP = pfit
err_BFP = np.array(error)

print ('Power Law: y = (%s+-%s)x  + (%s+-%s)' %(round(BFP[0],5),round(err_BFP[0],5),round(BFP[1],5),round(err_BFP[1],5)))

a1, a2 = BFP[0], BFP[1]
a1_err, a2_err = err_BFP[0], err_BFP[1]
zs = list(np.linspace(0.0051,20,10000))
scipy_fit = [rhoLya(z,a1,a2) for z in zs]
n = 1000
a1_list = np.random.normal(a1, a1_err*0.2, n)
a2_list = np.random.normal(a2, a2_err*0.2, n)

all_runs = []
for _ in a1_list:
    all_runs.append([])

for A,B,loop in zip(a1_list,a2_list,all_runs):
    for z in zs:
        parameter = np.random.randint(2)
        if parameter == 0:
            loop.append(rhoLya(z,A,a2))
        else:
            loop.append(rhoLya(z,a1,B))

median, median_upper_percentile, median_lower_percentile = [],[],[]
for z in zs:
    ind = zs.index(z)
    P_Lya = [i[ind] for i in all_runs]
    median_lower_percentile.append(np.percentile(P_Lya,16))
    median.append(np.median(P_Lya))
    median_upper_percentile.append(np.percentile(P_Lya,84))

#plt.plot(zs,scipy_fit,color='black')
plt.figure('LyaPowerLaw')
plt.fill_between(zs,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "grey", edgecolor = "black", linewidth = 1.2, label=r'$1\sigma$ conf. interval')
plt.plot(zs, median, label='ReHiLAE (this study, median)',color='black')
plt.scatter(xs,ys,color='black',marker='.',label='Sobral+2018')
plt.errorbar(xs,ys,y_sigma,color='black',ls='none')
plt.xlim(0,20)
plt.xlabel(r'Redshift ($z$)')
plt.ylabel(r'$log_{10}(\rho_{Ly\alpha}$ [$erg \ s^{-1} \ Mpc^{-3}$])')
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