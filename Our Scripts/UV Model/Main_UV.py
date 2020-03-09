import math
import numpy as np
import pylab
import matplotlib.pyplot as plt
from astropy.io import fits as py
from astropy.table import Table
from scipy.optimize import curve_fit
from scipy.integrate import odeint
from scipy.interpolate import interp1d


#CONSTANTS######################################################
H_0=0.0692 #Hubble's constant (Gyr⁻¹)
W_M=0.308 #Matter energy density parameter
W_lambda=0.692 #Dark energy density parameter
W_b_h_sqr=0.0223 #Baryon energy density parameter
X_p=0.75 #Hydrogen mass fraction
Y_p=0.25 #Helium mass fraction
f_esc_zero=2.3 #Escape fraction
Q_Hii_zero=0. #Q_Hii for Z=10
c_ha=1.36E-12 #Recombination coefficient
EW_avg=140.321828866 #Average equivalent width (angstroms)

alpha=1.17
T=20000 #K Temperature
C=3 #clumping factor
startT=0.124 #Start time
finishT=14 #Finish time
intervalNumber = 10000
TStep=(finishT - startT)/(intervalNumber) #Size of steps in time
EW = 148.9705
P1 = -5.288
P2 = 30.39894



def red(t,alt=False): #UNITLESS - calculates redshift from comsic time (Gyrs)
    if alt == False:
        return ((((28./(t))-1.)**(1./2.)-1.))
    else:
        return (((math.sinh(3*W_lambda**0.5*t / (2*1/H_0) ) * (W_lambda/W_M)**-0.5)**(-2/3)) )-1

def t(z, alt=False): #calculates comsic time (Gyrs) from redshift
    if alt == False:
        return (28.)/(z**2 + 2*z + 2)
    else:
        return ( 2 * math.pow(H_0, -1) ) / (3 * math.pow(W_lambda, 0.5) ) * math.sinh(  math.pow(W_lambda / W_M, 0.5) * math.pow(z + 1, -1.5) )

def redshift(startT, finishT, TStep, alt = False):
    
    T = np.arange(startT, finishT, TStep)
    Z=red(T)
    return Z,T

z,t = redshift(startT, finishT, TStep)



def alpha_beta(): #cm³ s⁻¹ - recombination coefficient
    return 2.6*(10**(-13)) * ((T/(10**4))**(-0.76))

def n_H(): #cm^-3  hydrogen number density
    return 1.67 * (10**(-7)) * (W_b_h_sqr / 0.02) * (X_p / 0.75)

def t_rec(z): #s recombination time
    return 0.23 #(alpha_beta() * n_H() * C * (1 + Y_p/(4*X_p)) * (1 + z)**(3))**(-1)

def E_ion(z):#Hz/erg
    return 10**(24.4 + math.log10(1 + z)) #10**25.6 

def f_esc_UV(z): #escape fraction of UV
        return (f_esc_zero*((1+z)/3)**alpha)/100 # UV escape fraction 0.23

def P_UV(z, P1, P2):
    return (P1*math.log10(z+1) + P2)

def n_ion_dot_UV(z, P1, P2): #s⁻¹ cm⁻³ 
    return f_esc_UV(z) * E_ion(z) * (10**(P1*math.log10(z+1) + P2)) / (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3

def Q_Hii_dot_UV(z,Q_Hii,P1, P2): #s⁻¹	def Q_Hii_dot(z,Q_Hii): #s⁻¹
    return (((n_ion_dot_UV(z, P1, P2)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16) # conversion from Gyr^-1 to s^-1	    return (((n_ion_dot(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16)  

def dQ_dt_UV(Q, t,P1, P2):
    
    dQ_dt = Q_Hii_dot_UV(red(t), Q,P1, P2)
    return dQ_dt

#GENERATE Q ARRAY

def main_UV(arguements):
    Q = odeint(dQ_dt_UV, Q_Hii_zero, t, args=(arguements))
    Q[Q>1.0] = 1.0 # 100% HII
    Q[Q<0.0] = 0.0 # 100% HI

    return Q
