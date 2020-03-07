from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math

"""
SC4K data shown in Sobral et al. 2018 table C3

available at: https://arxiv.org/pdf/1712.04451.pdf

returns P_L_Lya
"""

plt.figure('capped rho')

############################################# Power Law
def Function(x,a1,a2):
	return a1*x + a2

xs = pylab.array([2.2,2.5,2.8,3.0,3.2, 3.3, 3.7, 4.1, 4.6, 4.8, 5.1, 5.3, 5.8])# ,7,10,11,12,13,14,16])                               
y1 = pylab.array([0.52, 0.74, 0.77, 0.88, 0.84, 0.85, 1.01, 0.87, 1.19, 1.12, 1.27, 1.08, 1.10]) #,.7,0.5,0.3,0.10,0.05,0.001,0.001])  # data from SC4K Sobral 
ys = [math.log10(i*10**40) for i in y1]
xs_log = pylab.array([math.log10(i+1) for i in xs])
y_sigma_up = [0.05, 0.08, 0.1, 0.1, 0.9, 0.1, 0.2, 0.15, 0.4, 0.37, 0.48, 0.21, 0.21]#, 0.5 ,0.5,0.5,0.5,0.5, 0.5,0.5]
y_sigma_down =[0.05, 0.07, 0.09, 0.09, 0.09, 0.09, 0.16, 0.11, 0.3, 0.27, 0.32, 0.16, 0.16]#, 0.5,0.5,0.5,0.5,0.5, 0.5,0.5]
y_sigma = [(y_sigma_up[i]+y_sigma_down[i])/2 for i in range(len(y_sigma_down))]
p2 = pylab.polyfit(xs_log, ys, 1.0)
p1 = pylab.polyfit(xs, ys, 2.0)
p = pylab.poly1d(p2)
params = p.coefficients # initial estimation for parameters

# Proper fit: curvefit uses the function provided and the data with errors to minimise the residuals
pfit, pcov = curve_fit(Function, xs_log,ys,p0=params,sigma=y_sigma)

# Pfit will now contain the paramentes in the format [P_0,Phase0,A_0,Av_mag] (best fit parameters)
# Let's calculate the errors now from the covariance matrix:
error = [] # Empty error array
for i in range(len(pfit)): # Let's sample each individual parameter used
  try:
    error.append(np.absolute(pcov[i][i])**0.5)
  except:
    error.append(0.00)
BFP = pfit  # Best fit parameters
err_BFP = np.array(error) # Errors on best fit parameters

###### Print out the results for the best fit
print ('Power Law: y = (%s+-%s)x  + (%s+-%s)' %(round(BFP[0],5),round(err_BFP[0],5),round(BFP[1],5),round(err_BFP[1],5)))
# Show these on the plot as well

a1, a2, = BFP[0], BFP[1]
xdata = [i for i in range(14)]

def rhoLya(z):
    x = math.log10(1+z)
    if z < 5.8 :
      return 1.13869*x + 39.18853
    elif z > 5.8 :
      return  -5.288*x + 44.5388067

densities = [rhoLya(z) for z in xdata]

plt.figure(1)
plt.plot(xdata, densities, color = 'steelblue', label = ('Power Law: y = (%s+-%s)x  + (%s+-%s)' %(round(BFP[0],5),round(err_BFP[0],5),round(BFP[1],5),round(err_BFP[1],5))))


############################################## Quadratic
def Quad(x, a1, a2, a3):
  return a1*x**2 + a2*x + a3

p = pylab.poly1d(p1)
params = p.coefficients # initial estimation for parameters

# Proper fit: curvefit uses the function provided and the data with errors to minimise the residuals
pfit, pcov = curve_fit(Quad, xs,ys,p0=params,sigma=y_sigma)

# Pfit will now contain the paramentes in the format [P_0,Phase0,A_0,Av_mag] (best fit parameters)
# Let's calculate the errors now from the covariance matrix:
error = [] # Empty error array
for i in range(len(pfit)): # Let's sample each individual parameter used
  try:
    error.append(np.absolute(pcov[i][i])**0.5)
  except:
    error.append(0.00)
BFP = pfit  # Best fit parameters
err_BFP = np.array(error) # Errors on best fit parameters

###### Print out the results for the best fit
print ('Quadratic: y = (%s+-%s)$x^2$ + (%s+-%s)x + (%s+-%s)' %(round(BFP[0],5),round(err_BFP[0],5),round(BFP[1],5),round(err_BFP[1],5), round(BFP[2],5),round(err_BFP[2],5)))
# Show these on the plot as well

a1, a2, a3 = BFP[0], BFP[1], BFP[2], 
xdata = [i for i in range(0,14)]
data = [Quad(i,a1,a2, a3) for i in xdata]


plt.figure(1)
plt.scatter(xs, ys, color='black')
plt.errorbar(xs, ys, yerr=y_sigma, ls = 'none', color='black')
plt.plot(xdata, data, color='black', label=('Quadratic: y =(%s+-%s)$x^2$ + (%s+-%s)x + (%s+-%s)' %(round(BFP[0],3),round(err_BFP[0],3),round(BFP[1],3),round(err_BFP[1],3), round(BFP[2],3),round(err_BFP[2],3))))
plt.xlabel(r'Redshift (z)')
plt.ylabel(r'$log(\rho_{L_{Ly\alpha}}) [erg s^{−1} Mpc^{−3}]$')
plt.legend()
plt.show()


print(rhoLya(5.8))
