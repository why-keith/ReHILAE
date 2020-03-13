import matplotlib.pyplot as plt
import math
import numpy as np
import pylab
from scipy.optimize import curve_fit

def cubic(z,a1,a2,a3,a4):
    return a1*z**3 + a2*z**2 +a3*z +a4

def powerLaw(x,a5,a6):
    return a5*x + a6

def rho_Lya(z,a1,a2,a3,a4,a5,a6,scale):
    if z <=5.8:
        x = math.log10(1+z)
        return a5*x + a6
    if z>5.8:
        uv = a1*z**3 + a2*z**2 + a3*z + a4
        scaled = uv + scale
        return scaled

########################################## UV density
redshift_1 = pylab.array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4, 14])
UV_densities = pylab.array([26.52, 26.30, 26.10, 25.98, 25.67, 24.62, 23.00])
y_sigma = pylab.array([0.06,0.06,0.06,0.06,0.06, 0.04,0.001])                                 
cubic_fit = pylab.polyfit(redshift_1,UV_densities,3.0)
cubic_equation = pylab.poly1d(cubic_fit)
cubic_params = cubic_equation.coefficients

pfit, pcov = curve_fit(cubic, redshift_1, UV_densities, p0=cubic_params, sigma=y_sigma)
a1,a2,a3,a4 = pfit[0], pfit[1], pfit[2], pfit[3]

########################################## Lya density
x = [2.2,2.5,2.8,3.0,3.2, 3.3, 3.7, 4.1, 4.6, 4.8, 5.1, 5.3, 5.8 ]
xs = pylab.array([math.log10(i+1) for i in x])
y = [0.52, 0.74, 0.77, 0.88, 0.84, 0.85, 1.01, 0.87, 1.19, 1.12, 1.27, 1.08, 1.10] # data from SC4K Sobral 
ys = [math.log10(i*10**40) for i in y]
err = 0.15*np.array([0.52, 0.74, 0.77, 0.88, 0.84, 0.85, 1.01, 0.87, 1.19, 1.12, 1.27, 1.08, 1.10])
p2 = pylab.polyfit(xs, ys, 1.0)
p = pylab.poly1d(p2)
params = p.coefficients # initial estimation for parameters

# Proper fit: curvefit uses the function provided and the data with errors to minimise the residuals
pfit, pcov = curve_fit(powerLaw, xs,ys,p0=params,sigma=err)
a5,a6 = pfit[0], pfit[1]

########################################## Plotting
Lya_end = powerLaw(math.log10(1+5.8),a5,a6)
UV_start = cubic(5.8,a1,a2,a3,a4)
scale = Lya_end - UV_start


print('UV model: P_UV = %sz**3 + %sz**2 + %sz + %s' % (a1,a2,a3,a4))
print('Before z=5.8 : log(P_Lya) = %slog(1+z) + %s' % (a5,a6))
print('After z=5.8 : log(P_Lya) = %sz**3 + %sz**2 + %sz + %s + %s' %(a1,a2,a3,a4,scale))

zs = np.linspace(0.0051,12,100000)
UV = [cubic(z,a1,a2,a3,a4) for z in zs]
scaled_UV = [(i+scale) for i in UV]
Lya = [powerLaw(math.log10(1+z),a5,a6) for z in zs]
new_Lya = [rho_Lya(z,a1,a2,a3,a4,a5,a6,scale) for z in zs]

plt.figure()
plt.plot(zs,scaled_UV,color='blue',label='Scaled UV')
plt.plot(zs,Lya,color='red',label=r'Ly$\alpha$')
plt.plot(zs, new_Lya, color='green', linestyle='--',label=r'Combined Ly$\alpha$')
plt.xlabel('Redshift')
plt.ylabel(r'$log_{10}(\rho_{Ly\alpha}$ [$erg \ s^{-1} \ Mpc^{-3}$])')
plt.legend()

plt.figure()
plt.plot(zs,new_Lya,color='steelblue')
plt.xlabel('Redshift')
plt.ylabel(r'$log_{10}(\rho_{Ly\alpha}$ [$erg \ s^{-1} \ Mpc^{-3}$])')
plt.scatter(x,ys,color='black',marker='.')
plt.errorbar(x,ys, yerr=err, color='black', ls='none')

plt.show()