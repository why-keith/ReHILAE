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


def f_esc(z):
    return ((1+z)/3)**alpha * f_esc0

def P_UV(z):
    log_P_UV = C1*z**3 + C2*z**2 + C3*z + C4 
    return 10**(log_P_UV)

def nion(z):
    return f_esc(z)*xi_ion*P_UV(z)/ (2.938e+73)

def dQ_dt(Q,t):
    z= ((((28./(t))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift
    return ((nion(z))/n_H-Q/t_rec(z))*3.1536e+16 # conversion from Gyr^-1 to s^-1

def main(ts, arguements):
    global C1, C2, C3, C4,f1,f2 
    C1, C2, C3, C4, f1, f2 = arguements
    Q = odeint(dQ_dt, 0, ts,)
    Q[Q>1.0] = 1.0 # 100% HII
    Q[Q<0.0] = 0.0 # 100% HI
    return Q
