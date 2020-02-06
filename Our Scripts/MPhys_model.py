#Based on pg19 of MPhys Thesis
import math
import numpy as np
import pylab

#CONSTANTS
H_0=0.0692 #Gyr⁻¹ - Hubble's constant
W_M=0.308 #matter energy density parameter
W_lambda=0.692 #dark energy density parameter
W_b_h_sqr=0.0223 # baryon energy density parameter
X_p=0.75 #hydrogen mass fraction
Y_p=0.25 #helium mass fraction
f_esc_zero=2.3 #escape fraction zero
Q_Hii_zero=0.

#PARAMETER DEFAULTS
alpha=1.17
T=20000 #K Temperature
C=3 #clumping factor
startTime=0.051
finishTime=14
intervalNumber = 10000
timeStep = (finishTime - startTime)/(intervalNumber)
#########


#FUNCTIONS
def set_variables(_alpha=alpha,_T=T,_C=C): #allows changing of parameters from outside model.py
    global alpha
    global T
    global C
    
    alpha=_alpha
    T=_T
    C=_C
    
    return "α={} \nT={}\nC={}".format(alpha,T,C)
    
def z(t,alt=False): #UNITLESS - calculates redshift from comsic time (Gyrs)
    if alt == False:
        return ((((28./(t))-1.)**(1./2.)-1.))
    else:
        return (((math.sinh(3*W_lambda**0.5*t / (2*1/H_0) ) * (W_lambda/W_M)**-0.5)**(-2/3)) )-1

def t(z): #calculates comsic time (Gyrs) from redshift
    return ( 2 * math.pow(H_0, -1) ) / (3 * math.pow(W_lambda, 0.5) ) * math.sinh(  math.pow(W_lambda / W_M, 0.5) * math.pow(z + 1, -1.5) )

def alpha_beta(): #cm³ s⁻¹ - recombination coefficient
    return 2.6*(10**(-13)) * ((T/(10**4))**(-0.76))

def n_H(): #cm^-3  hydrogen number density
    return 1.67 * (10**(-7)) * (W_b_h_sqr / 0.02) * (X_p / 0.75)

def t_rec(z): #s recombination time
    return (alpha_beta() * n_H() * C * (1 + Y_p/(4*X_p)) * (1 + z)**(3))**(-1)

def f_esc(z): #escape fraction 
    return (f_esc_zero*((1+z)/3)**alpha)/100

def E_ion(z):#Hz/erg
    return 10**(24.4 + math.log10(1 + z))

def P_uv(z):   #Hz^-1 s^-1 Mpc^-3  UV luminosity density 
    x = pylab.array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4, 14])
    y = pylab.array([26.52, 26.30, 26.10, 25.98, 25.67, 24.62, 23.00])
    #0.0,0.45, 0.9, 1.3,1.8, 2.5 ,
    #10.0**25.72, 10**25.87, 10**26.05, 10**26.30, 10**26.32, 10**26.36,                                        
    p1 = pylab.polyfit(x,y,3.0)
    p = pylab.poly1d(p1)
    if z==14:
       return 10.0**26.52
    else:
       return 10**p(z)
#OLD
#def f_esc_Lya(EW_0):
 #   return 0.0048*EW_0 
#OLD
#def Q_ion_Lya(L_Lya, f_esc_Lya, EW_0):
 #   return L_Lya / (1 - f_esc_Lya)*(0.042 * EW_0)

def n_ion_dot(z): #s⁻¹ cm⁻³ 
    return f_esc(z) * E_ion(z) * P_uv(z) / (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3
#OLD
#def n_ion_dot_Lya(L_Lya, EW_0):
    return Q_ion_Lya * f_esc_Lya

def Q_Hii_dot(z,Q_Hii): #s⁻¹	def Q_Hii_dot(z,Q_Hii): #s⁻¹
    return (((n_ion_dot(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16) # conversion from Gyr^-1 to s^-1	    return (((n_ion_dot(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16)  
