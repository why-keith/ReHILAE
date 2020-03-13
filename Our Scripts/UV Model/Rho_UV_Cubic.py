from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math
# DATA

def Function(x,a1,a2,a3,a4):
    return a1*x**3. + a2*x**2. + a3*x + a4
"""
xs = pylab.array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4,17])
y = pylab.array([10.**(26.52), 10.**(26.30), 10.**(26.10), 10.**(25.98), 10.**(25.67), 10.**(24.62),10.**(23.12) ])
ys = pylab.array([math.log10(i) for i in y])
err = pylab.array([math.log10(0.1*i) for i in y])
"""
xs = pylab.array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4, 14])
ys = pylab.array([26.52, 26.30, 26.10, 25.98, 25.67, 24.62, 23.00])
err = pylab.array([0.06,0.06,0.06,0.06,0.06, 0.04,0.001]) 
p1 = pylab.polyfit(xs,ys,3.)
p = pylab.poly1d(p1) # rho UV
params = p.coefficients # initial estimation for parameters

# Proper fit: curvefit uses the function provided and the data with errors to minimise the residuals
pfit, pcov = curve_fit(Function, xs,ys,p0=params,sigma=err)

error = [] # Empty error array
for i in range(len(pfit)): # Let's sample each individual parameter used
  try:
    error.append(np.absolute(pcov[i][i])**0.5)
  except:
    error.append(0.00)
BFP = pfit  # Best fit parameters
err_BFP = np.array(error) # Errors on best fit parameters

print ('y = (%s+-%s)x**3 + (%s+-%s)x**2 + (%s+-%s)x* + (%s+-%s)' %(round(BFP[0],7),round(err_BFP[0],7),round(BFP[1],7),round(err_BFP[1],7),round(BFP[2],7),round(err_BFP[2],7),round(BFP[3],7),round(err_BFP[3],7)))

a1, a2, a3, a4, = BFP[0], BFP[1], BFP[2], BFP[3]

zs = np.linspace(0.0051,14,100000)

a1_min, a2_min, a3_min, a4_min = BFP[0]-err_BFP[0], BFP[1]-err_BFP[1], BFP[2]-err_BFP[2], BFP[3]-err_BFP[3]
a1_max, a2_max, a3_max, a4_max = BFP[0]+err_BFP[0], BFP[1]+err_BFP[1], BFP[2]+err_BFP[2], BFP[3]+err_BFP[3] 
data = [Function(i,a1,a2,a3,a4) for i in zs]
data_min = [Function(i,a1_min,a2_min,a3_min,a4_min) for i in zs]
data_max = [Function(i,a1_max,a2_max,a3_max,a4_max) for i in zs]
polyfit_data = [p(i) for i in zs]

plt.plot(zs, data, color='red', label='Average optimsed data')
plt.plot(zs, polyfit_data, color='black', linestyle='-.', label='Polyfit data')
plt.plot(zs, data_min, color='blue', linestyle='--', label='Minimum optimsed data')
plt.plot(zs, data_max, color='green', linestyle='--',  label='Maximum optimsed data')
plt.xlabel('Redshift ($z$)')
plt.ylabel(r'$\log (\rho_{UV}) \ [erg \ s^{-1} \ Mpc^{-3}]$')
plt.title('Plot to show Scipy and polyfit curves')
plt.legend()
plt.show()