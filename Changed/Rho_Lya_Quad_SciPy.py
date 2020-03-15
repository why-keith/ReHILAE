from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math
# DATA
##########################################################################################################
def Function(x,a1,a2,a3):
	return a1*x**2.+a2*x+a3

xs = [2.2,2.5,2.8,3.0,3.2, 3.3, 3.7, 4.1, 4.6, 4.8, 5.1, 5.3, 5.8 ]
y = np.array([0.52, 0.74, 0.77, 0.88, 0.84, 0.85, 1.01, 0.87, 1.19, 1.12, 1.27, 1.08, 1.10]) # data from SC4K Sobral 
ys = np.array([math.log10(i*10**40) for i in y])
ys_err = 0.15*y
p2 = pylab.polyfit(xs, ys, 2.0)
p = pylab.poly1d(p2)
params = p.coefficients
print(p) # initial estimation for parameters

# Proper fit: curvefit uses the function provided and the data with errors to minimise the residuals
pfit, pcov = curve_fit(Function, xs,ys,p0=params,sigma=ys_err)

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
print (r'$\log(\rho_{L_{Ly\alpha}})$ =(%s+-%s)$z^{2}$ + (%s+-%s)$z$+ (%s+-%s)' %(round(BFP[0],2),round(err_BFP[0],2),round(BFP[1],2),round(err_BFP[1],2),round(BFP[2],2),round(err_BFP[2],2)))
# Show these on the plot as well

a1, a2, a3, = BFP[0], BFP[1], BFP[2]
xdata = [i for i in range(0,14)]
data = [Function(i,a1,a2,a3) for i in range(0,14)]
data_min = [Function(i,a1-err_BFP[0],a2-err_BFP[1],a3-err_BFP[2]) for i in xdata]
data_max = [Function(i,a1+err_BFP[0],a2+err_BFP[1],a3+err_BFP[2]) for i in xdata]

################################################################################
plt.figure("Rho_Lya_Quad")
plt.scatter(xs, ys, color='black', marker='.')
plt.errorbar(xs, ys, yerr=ys_err, ls = 'none', color='black')
plt.plot(xdata, data, color='steelblue')
plt.fill_between(xdata, data_min, data_max, alpha=0.4, color = "grey", edgecolor = "black", linewidth = 1.2)
plt.ylabel(r'$log_{10}(\rho_{Ly\alpha} \ [erg \ s^{-1} \ Mpc^{-3})]$')
plt.xlabel(r'Redshift ($z$)')
plt.show()