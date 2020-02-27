import math
import numpy as np
import pylab
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import curve_fit

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



#PARAMETER DEFAULTS#############################################
alpha=1.17
T=20000 #K Temperature
C=3 #clumping factor
startT=0.124 #Start time
finishT=14 #Finish time
intervalNumber = 10000
TStep=(finishT - startT)/(intervalNumber) #Size of steps in time
EW = 148.9705

def alpha_beta(): #cm³ s⁻¹ - recombination coefficient
    return 2.6*(10**(-13)) * ((T/(10**4))**(-0.76))

def n_H(): #cm^-3  hydrogen number density
    return 1.67 * (10**(-7)) * (W_b_h_sqr / 0.02) * (X_p / 0.75)

def t_rec(z): #s recombination time
    return (alpha_beta() * n_H() * C * (1 + Y_p/(4*X_p)) * (1 + z)**(3))**(-1)

def P_L_Lya(x, P1, P2, P3):
    return P1*x**2 + P2*x +P3

def f_esc_LyC(x, f1, f2):
    return f1*x + f2

def Q_ion_LyC(z, P1, P2, P3, f1, f2): # replaces P_uv and E_ion	
    return P_L_Lya(z, P1, P2, P3) / ((c_ha*(1 - f_esc_LyC(EW, f1, f2), f1, f2))*(0.042 * EW))

def n_ion_dot_LyC(z, P1, P2, P3,  f1, f2): # replaces n_ion_dot using Q_ion_LyC	
    if Q_ion_LyC(z, P1, P2, P3,  f1, f2) * f_esc_LyC(EW, f1, f2) / (2.938e+73) <= 0 :
        return 0 
    else:
        return Q_ion_LyC(z, P1, P2, P3, f1, f2) * f_esc_LyC(EW, f1, f2) / (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3

def Q_Hii_dot(z, P1, P2, P3, f1, f2,Q_Hii): #s⁻¹	def Q_Hii_dot(z,Q_Hii): #s⁻¹
    return (((n_ion_dot_LyC(z, P1, P2, P3, f1, f2)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16) # conversion from Gyr^-1 to s^-1	    return (((n_ion_dot(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16)  
