from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math
import random_array_generator as rag

###############################################################################
xs = pylab.array([2.2,2.5,2.8,3.0,3.2, 3.3, 3.7, 4.1, 4.6, 4.8, 5.1, 5.3, 5.8])# ,7,10,11,12,13,14,16])                               
y1 = pylab.array([0.52, 0.74, 0.77, 0.88, 0.84, 0.85, 1.01, 0.87, 1.19, 1.12, 1.27, 1.08, 1.10]) #,.7,0.5,0.3,0.10,0.05,0.001,0.001])  # data from SC4K Sobral 
ys = [math.log10(i*10**40) for i in y1]
xs_log = pylab.array([math.log10(i+1) for i in xs])
y_sigma_up = [0.05, 0.08, 0.1, 0.1, 0.9, 0.1, 0.2, 0.15, 0.4, 0.37, 0.48, 0.21, 0.21]#, 0.5 ,0.5,0.5,0.5,0.5, 0.5,0.5]
y_sigma_down =[0.05, 0.07, 0.09, 0.09, 0.09, 0.09, 0.16, 0.11, 0.3, 0.27, 0.32, 0.16, 0.16]#, 0.5,0.5,0.5,0.5,0.5, 0.5,0.5]
y_sigma = [(y_sigma_up[i]+y_sigma_down[i])/2 for i in range(len(y_sigma_down))]

def rhoLya(z):
    x = math.log10(1+z)
    if z < 5.8 :
      return 1.13869*x + 39.18853
    elif z > 5.8 :
      return  -5.288*x + 44.5388067 

xdata = [i for i in range(0, 14)]
densities = [rhoLya(z) for z in xdata]
###############################################################################
P1_1 = [1.13869]
P1_1_error = [0.18238]

P1_2 = [39.18853]
P1_2_error = [0.10709]

P2_1 = [-5.288]
P2_1_error = [0.83404]

P2_2 = [44.5388067]
P2_2_error = [0.757]


P1_1_P_L_Lya = rag.random_Arrays(len(P1_1),P1_1,P1_1_error,P1_1_error)
P1_2_P_L_Lya = rag.random_Arrays(len(P1_2),P1_2 ,P1_2_error,P1_2_error)
P2_1_P_L_Lya = rag.random_Arrays(len(P2_1),P2_1,P2_1_error,P2_1_error)
P2_2_P_L_Lya = rag.random_Arrays(len(P2_2),P2_2 ,P2_2_error,P2_2_error)

p_lya_iterations = []
for i in range(len(P1_1_P_L_Lya)):
    dummy = []
    for j in xdata:
        if j <5.8:
            dummy.append(P1_1_P_L_Lya[i][0]*math.log10(1+j) + P1_2_P_L_Lya[i][0])
        else:
            dummy.append(P2_1_P_L_Lya[i][0]*math.log10(1+j) + P2_2_P_L_Lya[i][0])
    p_lya_iterations.append(dummy)
        
median, median_lower_percentile, median_upper_percentile = rag.median_y_values(len(p_lya_iterations[0]),p_lya_iterations)


###############################################################################
plt.figure("Rho_corr")
plt.scatter(xs, ys, color='black')
plt.xlabel(r'Redshift (z)')
plt.ylabel(r'$log(\rho_{L_{Ly\alpha}}) [erg s^{−1} Mpc^{−3}]$')
plt.errorbar(xs, ys, yerr=y_sigma, ls = 'none', color='black')
plt.plot(xdata, densities, color = 'steelblue')

plt.plot(xdata, median, "--")
plt.fill_between(xdata,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "grey", edgecolor = "black", linewidth = 1.2)

plt.show()