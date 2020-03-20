from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math
import matplotlib

def Function(x,a1,a2,a3,a4):
    return a1*x**3. + a2*x**2. + a3*x + a4

xs = pylab.array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4, 14])
ys = pylab.array([26.52, 26.30, 26.10, 25.98, 25.67, 24.62, 23.00])
err = pylab.array([0.06,0.06,0.06,0.06,0.06, 0.04,0.001]) 
p1 = pylab.polyfit(xs,ys,3.)
p = pylab.poly1d(p1) 
params = p.coefficients

pfit, pcov = curve_fit(Function,xs,ys,p0=params,sigma=err)

error = [] 
for i in range(len(pfit)): 
  try:
    error.append(np.absolute(pcov[i][i])**0.5)
  except:
    error.append(0.00)
BFP = pfit  
err_BFP = np.array(error)

print ('log(P_UV) = (%s+-%s)z**3 + (%s+-%s)z**2 + (%s+-%s)x + (%s+-%s)' %(round(BFP[0],7),round(err_BFP[0],7),round(BFP[1],7),round(err_BFP[1],7),round(BFP[2],7),round(err_BFP[2],7),round(BFP[3],7),round(err_BFP[3],7)))

a1, a2, a3, a4, = BFP[0], BFP[1], BFP[2], BFP[3]
a1_err, a2_err, a3_err, a4_err = err_BFP[0], err_BFP[1], err_BFP[2], err_BFP[3]

zs = list(np.linspace(0.0051,20,10000))
scipy_fit = [Function(z,a1,a2,a3,a4) for z in zs]
number_of_iterations = 1000

a1_list = np.random.normal(a1, a1_err*0.2, number_of_iterations)
a2_list = np.random.normal(a2, a2_err*0.2, number_of_iterations)
a3_list = np.random.normal(a3, a3_err*0.2, number_of_iterations)
a4_list = np.random.normal(a4, a4_err*0.2, number_of_iterations)

all_runs = []
for _ in a1_list:
    all_runs.append([])

for A,B,C,D,loop in zip(a1_list,a2_list,a3_list,a4_list,all_runs):
    for z in zs:
        parameter = np.random.randint(4)
        if parameter == 0:
            loop.append(Function(z,A,a2,a3,a4))
        elif parameter == 1:
            loop.append(Function(z,a1,B,a3,a4))
        elif parameter == 2:
            loop.append(Function(z,a1,a2,C,a4))
        else:
            loop.append(Function(z,a1,a2,a3,D))

median, median_upper_percentile, median_lower_percentile = [],[],[]
for z in zs:
    ind = zs.index(z)
    LumDensities = [i[ind] for i in all_runs]
    median_lower_percentile.append(np.percentile(LumDensities,16))
    median.append(np.median(LumDensities))
    median_upper_percentile.append(np.percentile(LumDensities,84))

plt.figure("Rho_UV_Cubic")
#plt.plot(zs, scipy_fit, color='black', label='SciPy fit')
plt.errorbar(xs,ys,yerr=err,ls='none',color='black')
plt.fill_between(zs,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "grey", edgecolor = "black", linewidth = 1.2, label=r'$1\sigma$ conf. interval')
plt.plot(zs, median, label='ReHiLAE (this study, median)',color='black')
plt.scatter(xs,ys,color='black',marker='.',label='Bouwens+2015')
plt.xlabel(r'Redshift (z)')
plt.ylabel(r'$log_{10}(\rho_{UV} \ [erg \ s^{-1} \ Mpc^{-3}])$')
plt.xlim(0,20)
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