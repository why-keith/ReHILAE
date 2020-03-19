from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math
import matplotlib

a1,a2,a3,a4,a5,a6 = 1.13869,39.18853,0.0013679,-0.0559273,0.2937793,26.0952716
a1_err,a2_err,a3_err,a4_err,a5_err,a6_err = 0.18238,0.10709,0.0009035,0.0248473,0.2101038,0.5405448

redshift_1 = [2.2,2.5,2.8,3.0,3.2, 3.3, 3.7, 4.1, 4.6, 4.8, 5.1, 5.3, 5.8 ]
y = [0.52, 0.74, 0.77, 0.88, 0.84, 0.85, 1.01, 0.87, 1.19, 1.12, 1.27, 1.08, 1.10]
LyaDensities = [math.log10(i*10**40) for i in y]
Lya_err = [0.05,0.075,0.095,0.095,0.09,0.095,0.18,0.13,0.35,0.32,0.4,0.185,0.185]

redshift_2 = pylab.array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4, 14])
UV_densities = pylab.array([26.52, 26.30, 26.10, 25.98, 25.67, 24.62, 23.00])
UV_err = pylab.array([0.06,0.06,0.06,0.06,0.06, 0.04,0.001])

def RhoLyaPower(z,a1=a1,a2=a2):
    x = math.log10(1+z)
    return a1*x + a2

def RhoUV(z,a3=a3,a4=a4,a5=a5,a6=a6):
    return a3*z**3 + a4*z**2 + a5*z + a6

def rhoLya(z,a1=a1,a2=a2,a3=a3,a4=a4,a5=a5,a6=a6):
    if z <=5.8:
        x = math.log10(1+z)
        return a1*x + a2
    if z>5.8:
        uv = a3*z**3 + a4*z**2 + a5*z + a6
        scaled = uv + scale
        return scaled

UV_start = RhoUV(5.8)
Lya_end = RhoLyaPower(5.8)
scale = Lya_end - UV_start # 13.951808701009465

zs = list(np.linspace(0,20,10000))

scaledUV_min, scaledUV, scaledUV_max = [],[],[]
Lya_min, Lya, Lya_max = [],[],[]
for z in zs:
    scaledUV.append(scale+ RhoUV(z))
    Lya.append(RhoLyaPower(z))

scipy_fit = [rhoLya(z) for z in zs]
n = 1000
a1_list = np.random.normal(a1, a1_err*0.2, n)
a2_list = np.random.normal(a2, a2_err*0.2, n)
a3_list = np.random.normal(a3, a3_err*0.2, n)
a4_list = np.random.normal(a4, a4_err*0.2, n)
a5_list = np.random.normal(a5, a5_err*0.2, n)
a6_list = np.random.normal(a6, a6_err*0.2, n)

all_runs = []
for _ in a1_list:
    all_runs.append([])

for A,B,C,D,E,F,loop in zip(a1_list,a2_list,a3_list,a4_list,a5_list,a6_list,all_runs):
    for z in zs:
        if z <= 5.8:
            parameter = np.random.randint(2)
            if parameter == 0:
                loop.append(rhoLya(z,A,a2))
            else:
                loop.append(rhoLya(z,a1,B))
        else:
            parameter = np.random.randint(4)
            if parameter == 0:
                loop.append(rhoLya(z,a3=C))
            elif parameter == 1:
                loop.append(rhoLya(z,a4=D))
            elif parameter == 2:
                loop.append(rhoLya(z,a5=E))
            else:
                loop.append(rhoLya(z,a6=F))

median, median_upper_percentile, median_lower_percentile = [],[],[]
for z in zs:
    ind = zs.index(z)
    LumDensities = [i[ind] for i in all_runs]
    median_lower_percentile.append(np.percentile(LumDensities,16))
    median.append(np.median(LumDensities))
    median_upper_percentile.append(np.percentile(LumDensities,84))

plt.figure('Luminosity Densities')
plt.plot(zs,scaledUV,color='grey',label=r'Scaled $\rho_{UV}$')
plt.plot(zs,Lya,color='black', label=r'$\rho_{Ly\alpha}$')
plt.scatter(redshift_1,LyaDensities,color='black',marker='.',label='Sobral+2018')
plt.errorbar(redshift_1,LyaDensities, yerr=Lya_err, color='black', ls='none')
plt.scatter(redshift_2,[i+scale for i in UV_densities],color='blue',marker='.',label='Bouwens+2015')
plt.errorbar(redshift_2,[i+scale for i in UV_densities], yerr=UV_err, color='blue', ls='none')
plt.xlabel(r'Redshift ($z$)')
plt.ylabel(r'Luminosity Density [$erg \ s^{-1} \ Mpc^{-3}$]')
plt.xlim(0,20)
plt.tick_params(which='both',direction='in',right=True,top=True)
plt.legend()

plt.figure('Rho Lya Combined')
#plt.plot(zs,scipy_fit, color='black',label=r'SciPy fit')
plt.fill_between(zs,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "grey", edgecolor = "black", linewidth = 1.2, label=r'$1\sigma$ conf. interval')
plt.plot(zs, median, label='ReHiLAE (this study, median)',color='black')
plt.scatter(redshift_1,LyaDensities,color='black',marker='.',label='Sobral+2018')
plt.errorbar(redshift_1,LyaDensities, yerr=Lya_err, color='black', ls='none')
plt.scatter(redshift_2[4:],[i+scale for i in UV_densities[4:]],color='blue',marker='.',label='Bouwens+2015')
plt.errorbar(redshift_2[4:],[i+scale for i in UV_densities[4:]], yerr=UV_err[4:], color='blue', ls='none')
plt.xlabel(r'Redshift ($z$)')
plt.ylabel(r'$log_{10}(\rho_{Ly\alpha}$ [$erg \ s^{-1} \ Mpc^{-3}$])')
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