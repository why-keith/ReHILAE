from pylab import array, polyfit, poly1d
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

xi_ion = 10.**(25.6) #production_efficiency [Hz erg**(-1.)]
f_esc0 = 0.23
alpha = 1.17 # f_esc dependence on redshift
x1 = array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4,17])
y1 = array([10.**(26.52), 10.**(26.30), 10.**(26.10), 10.**(25.98), 10.**(25.67), 10.**(24.62),10.**(23.12) ])
p1 = polyfit(x1,y1,3.)
p = poly1d(p1) # rho UV
Omega_b = 0.04 #baryon_density
h = 0.7
X_p = 0.75 #hydrogen_mass_fraction
T = 20000. #IGM_temperature [K]
C = 3  #clumping_factor
Y_p = (1.-X_p) #helium mass fraction

alpha_B = 2.6*pow(10,-13.)*(T/10.**(4.))**(-0.76)  #case_B_recombination_coefficient
n_H = 1.67*pow(10,-7.) * ((Omega_b * h* h)/0.02) * (X_p/0.75) # per cm^3

def dQ_dt(Q,t):
    z= ((((28./(t))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift
    
    trec = 1./(C*alpha_B*(1.+Y_p/(4.*X_p))*n_H*(1.+z)**(3.)) # units: seconds 
    nion = f_esc0*xi_ion*p(z)/ (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3
    return ((nion)/n_H-Q/trec)*3.1536e+16 # conversion from Gyr^-1 to s^-1

ts = np.linspace(0.051,14,10000000) # time in Gyr
zs= ((((28./(ts))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift

Q = odeint(dQ_dt, 0., ts) 
Q[Q>1.0] = 1.0 # 100% HII boundary condition up
Q[Q<0.0] = 0.0 # 100% HI boundary condition down

plt.plot(zs, Q,)
plt.title('Reproduced Model')
plt.xlabel('Redshift ($z$)')
plt.ylabel('Fraction of Ionised Hydrogen ($Q_{II}$)')
plt.show()