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
EW = 148.9705
c_ha=1.36E-12

alpha_B = 2.6*pow(10,-13.)*(T/10.**(4.))**(-0.76)  #case_B_recombination_coefficient
n_H = 1.67*pow(10,-7.) * ((Omega_b * h* h)/0.02) * (X_p/0.75) # per cm^3


init_conditions={"C1":np.array([]), "C2":np.array([]), "C3":np.array([]), "C4":np.array([]), "P1":np.array([]),"P2":np.array([]),"F1":np.array([]),"F2":np.array([])} #records the inital conditions of each iteration

def t_rec(z):
    return 1./(C*alpha_B*(1.+Y_p/(4.*X_p))*n_H*(1.+z)**(3.)) # units: seconds

def f_esc():
    return F1*EW + F2

def P_Lya(z):
    if z <=5.8:
        return P1*math.log10(1+z) + P2
    if z>5.8:
        log_P_UV = C1*z**3 + C2*z**2 + C3*z + C4
        log_scaled = log_P_UV + 13.933143098531438
        return log_scaled

def Q_dot(z):
    return 10**P_Lya(z) / (c_ha*(1-f_esc())*(0.0042*EW))

def nion(z):
    return f_esc()*Q_dot(z)/ (2.938e+73)

def dQ_dt(Q,t):
    z= ((((28./(t))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift
    return ((nion(z))/n_H-Q/t_rec(z))*3.1536e+16 # conversion from Gyr^-1 to s^-1

def main(ts, arguements):
    global C1, C2, C3, C4, P1, P2, F1, F2
    C1, C2, C3, C4, P1, P2, F1, F2 = arguements
    Q = odeint(dQ_dt, 0, ts,)
    Q[Q>1.0] = 1.0 # 100% HII
    Q[Q<0.0] = 0.0 # 100% HI
    for i,j in zip(init_conditions,range(len(arguements))):        
        init_conditions[i]=np.append(init_conditions[i],arguements[j])
    
    return Q
