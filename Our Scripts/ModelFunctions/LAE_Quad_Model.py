import numpy as np
from scipy.integrate import odeint
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
c_ha=1.36E-12
scale = 13.951808701009465

alpha_B = 2.6*pow(10,-13.)*(T/10.**(4.))**(-0.76)  #case_B_recombination_coefficient
n_H = 1.67*pow(10,-7.) * ((Omega_b * h* h)/0.02) * (X_p/0.75) # per cm^3

def t_rec(z):
    return 1./(C*alpha_B*(1.+Y_p/(4.*X_p))*n_H*(1.+z)**(3.)) # units: seconds

def EW(z):
    return e1*z + e2

def f_esc(z):
    return f1*EW(z) + f2

def P_Lya(z):
    return p1*z**2 + p2*z + p3

def Q_dot(z):
    return 10**P_Lya(z) / (c_ha*(1-f_esc(z))*(0.0042*EW(z)))

def nion(z):
    return f_esc(z)*Q_dot(z)/ (2.938e+73)

def dQ_dt(Q,t):
    z= ((((28./(t))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift
    return ((nion(z))/n_H-Q/t_rec(z))*3.1536e+16 # conversion from Gyr^-1 to s^-1

def main(ts, arguements):
    global e1, e2, f1, f2, p1, p2, p3
    e1, e2, f1, f2, p1, p2, p3 = arguements
    Q = odeint(dQ_dt, 0, ts,)
    Q[Q>1.0] = 1.0 # 100% HII
    Q[Q<0.0] = 0.0 # 100% HI
    return Q