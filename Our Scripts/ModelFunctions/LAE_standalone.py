from pylab import array, polyfit, poly1d
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math

xi_ion = 10**25.4
f_esc0 = 0.23
alpha = 1.17 # f_esc dependence on redshift
Omega_b = 0.04 #baryon_density
h = 0.7
X_p = 0.75 #hydrogen_mass_fraction
T = 20000. #IGM_temperature [K]
C = 3  #clumping_factor
Y_p = (1.-X_p) #helium mass fraction
EW = 148.9705
c_ha=1.36E-12

alpha_B = 2.6*pow(10,-13.)*(T/10.**(4.))**(-0.76)  #case_B_recombination_coefficient
n_H = 1.67*pow(10,-7.) * ((Omega_b * h* h)/0.02) * (X_p/0.75) # per cm^3

def t_rec(z):
    return 1./(C*alpha_B*(1.+Y_p/(4.*X_p))*n_H*(1.+z)**(3.)) # units: seconds

def f_esc(z):
    return 0.1

def P_Lya(z):
    if z <=5.8:
        return 1.0388337446404217*math.log10(1+z) + 39.2529986910512
    if z>5.8:
        log_P_UV = 0.0013679194302549992*z**3 + -0.05592731916980156*z**2 + 0.29377928649363305*z + 26.095271603086044
        log_scaled = log_P_UV + 13.933143098531438
        return log_scaled

def Q_dot(z):
    return 10**P_Lya(z) / (c_ha*(1-f_esc(z))*(0.0042*EW))

def nion(z):
    return f_esc(z)*Q_dot(z)/ (2.938e+73)

def dQ_dt(Q,t):
    z= ((((28./(t))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift
    return ((nion(z))/n_H-Q/t_rec(z))*3.1536e+16 # conversion from Gyr^-1 to s^-1

ts = np.linspace(0.051,14,10000000) # time in Gyr
zs= ((((28./(ts))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift

Q = odeint(dQ_dt, 0., ts) 
Q[Q>1.0] = 1.0 # 100% HII boundary condition up
Q[Q<0.0] = 0.0 # 100% HI boundary condition down

plt.plot(zs, Q,)
plt.title('LAE Model')
plt.xlabel('Redshift ($z$)')
plt.ylabel('Fraction of Ionised Hydrogen ($Q_{II}$)')
plt.show()