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

error_1 = [] # Empty error array
for i in range(len(pfit)): # Let's sample each individual parameter used
  try:
    error_1.append(np.absolute(pcov[i][i])**0.5)
  except:
    error_1.append(0.00)
BFP_1 = pfit  # Best fit parameters
err_BFP_1 = np.array(error_1) # Errors on best fit parameters

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

error_2 = [] # Empty error array
for i in range(len(pfit)): # Let's sample each individual parameter used
  try:
    error_2.append(np.absolute(pcov[i][i])**0.5)
  except:
    error_2.append(0.00)
BFP = pfit  # Best fit parameters
err_BFP_2 = np.array(error_2) # Errors on best fit parameters

########################################## Plotting
Lya_end = powerLaw(math.log10(1+5.8),a5,a6)
UV_start = cubic(5.8,a1,a2,a3,a4)
scale = Lya_end - UV_start


print('UV model: log(P_UV) = (%s+%s)z**3 + (%s±%s)z**2 + (%s±%s)z + (%s±%s)' % (
    round(a1,5), round(err_BFP_1[0],5),
    round(a2,5), round(err_BFP_1[1],5),
    round(a3,5), round(err_BFP_1[2],5),
    round(a4,5), round(err_BFP_1[3],5)))
print('Before z=5.8 : log(P_Lya) = (%s±%s)log(1+z) + (%s±%s)' % (round(a5,5), round(err_BFP_2[0],5),round(a5,5), round(err_BFP_1[1],5)))
print('After z=5.8 : log(P_Lya) = (%s+%s)z**3 + (%s±%s)z**2 + (%s±%s)z + (%s±%s) + %s' %(
    round(a1,5), round(err_BFP_1[0],5),
    round(a2,5), round(err_BFP_1[1],5),
    round(a3,5), round(err_BFP_1[2],5),
    round(a4,5), round(err_BFP_1[3],5),
    round(scale,5)))

zs = np.linspace(0.0051,14,100000)

UV = [cubic(z,a1,a2,a3,a4) for z in zs]
UV_min = [cubic(z,a1-err_BFP_1[0],a2-err_BFP_1[1],a3-err_BFP_1[2],a4-err_BFP_1[3]) for z in zs]
UV_max = [cubic(z,a1+err_BFP_1[0],a2+err_BFP_1[1],a3+err_BFP_1[2],a4+err_BFP_1[3]) for z in zs]
scaled_UV = [(i+scale) for i in UV]
scaled_UV_min = [(i+scale) for i in UV_min]
scaled_UV_max = [(i+scale) for i in UV_max]

Lya = [powerLaw(math.log10(1+z),a5,a6) for z in zs]
Lya_min = [powerLaw(math.log10(1+z),a5-err_BFP_2[0],a6-err_BFP_2[1]) for z in zs]
Lya_max = [powerLaw(math.log10(1+z),a5+err_BFP_2[0],a6+err_BFP_2[1]) for z in zs]


plt.figure()
plt.plot(zs,scaled_UV,color='blue',label='Scaled UV density')
plt.plot(zs,Lya,color='red',label=r'Ly$\alpha$')
plt.scatter(redshift_1, [i+scale for i in UV_densities], color='black', marker='.')
plt.errorbar(redshift_1, [i+scale for i in UV_densities], yerr=y_sigma, ls='none', color='black')
plt.scatter(x,ys,color='black',marker='.')
plt.errorbar(x,ys, yerr=err, color='black', ls='none')
#plt.fill_between(zs, scaled_UV_min, scaled_UV_max, color='steelblue', alpha=0.5)
#plt.fill_between(zs, Lya_min, Lya_max, color='salmon', alpha=0.5)
plt.xlabel('Redshift')
plt.ylabel(r'$log_{10}(\rho_{Ly\alpha}$ [$erg \ s^{-1} \ Mpc^{-3}$])')
plt.legend()

new_Lya = [rho_Lya(z,a1,a2,a3,a4,a5,a6,scale) for z in zs]

plt.figure()
plt.plot(zs,new_Lya,color='steelblue')
plt.xlabel('Redshift')
plt.ylabel(r'$log_{10}(\rho_{Ly\alpha}$ [$erg \ s^{-1} \ Mpc^{-3}$])')
plt.scatter(x,ys,color='black',marker='.')
plt.errorbar(x,ys, yerr=err, color='black', ls='none')
plt.scatter(redshift_1[4:], [i+scale for i in UV_densities[4:]], color='black', marker='.')
plt.errorbar(redshift_1[4:], [i+scale for i in UV_densities[4:]], yerr=y_sigma[4:], ls='none', color='black')

plt.show()