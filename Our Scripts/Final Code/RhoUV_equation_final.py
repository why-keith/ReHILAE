from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math
# DATA

def Function(x,a1,a2):
	return a1*x + a2

xs = pylab.array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4])
ys = pylab.array([26.52, 26.30, 26.10, 25.98, 25.67, 24.62])
xs_log = pylab.array([math.log10(i+1) for i in xs])
y_sigma = pylab.array([0.06,0.06,0.06,0.06,0.06, 0.04])   

p2 = pylab.polyfit(xs_log, ys, 1.0)
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
print ('Function: y = (%s+-%s)x  + (%s+-%s)' %(round(BFP[0],5),round(err_BFP[0],5),round(BFP[1],5),round(err_BFP[1],5)))
# Show these on the plot as well

a1, a2, = BFP[0], BFP[1]
xdata = [i for i in range(0,round(math.log(15)))]
data = [Function(i,a1,a2) for i in xdata]

plt.figure(1)
plt.scatter(xs_log, ys, color='black')
#plt.errorbar(xs, ys, yerr=y_sigma, ls = 'none', color='black')
plt.plot(xdata, data, color='red', label=r'foo')
plt.xlabel(r'$log(z+1)$')
plt.ylabel(r'$log(\rho_{UV})$')
plt.legend()
plt.show()

def rhoUV(z):
    z = math.log10(1+z)
    return -5.288*z + 30.39894

zs = np.linspace(0.051, 14,10000)
densities = [rhoUV(z) for z in zs]

plt.plot(zs, densities)
plt.xlabel('Redshift (z)')
plt.ylabel(r'$log(\rho_{UV})$')
plt.show()

print(rhoUV(5.8))