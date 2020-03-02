"""
Defines, calculates and allows for the changing of variables used in the toy model 
"""
import math
import numpy as np
import pylab
from astropy.io import fits
from astropy.table import Table



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

def z(t,alt=False): #UNITLESS - calculates redshift from comsic time (Gyrs)
    if alt == False:
        return ((((28./(t))-1.)**(1./2.)-1.))
    else:
        return (((math.sinh(3*W_lambda**0.5*t / (2*1/H_0) ) * (W_lambda/W_M)**-0.5)**(-2/3)) )-1

def t(z, alt=False): #calculates comsic time (Gyrs) from redshift
    if alt == False:
        return (28.)/(z**2 + 2*z + 2)
    else:
        return ( 2 * math.pow(H_0, -1) ) / (3 * math.pow(W_lambda, 0.5) ) * math.sinh(  math.pow(W_lambda / W_M, 0.5) * math.pow(z + 1, -1.5) )

def alpha_beta(): #cm³ s⁻¹ - recombination coefficient
    return 2.6*(10**(-13)) * ((T/(10**4))**(-0.76))

def n_H(): #cm^-3  hydrogen number density
    return 1.67 * (10**(-7)) * (W_b_h_sqr / 0.02) * (X_p / 0.75)

def t_rec(z): #s recombination time
    return (alpha_beta() * n_H() * C * (1 + Y_p/(4*X_p)) * (1 + z)**(3))**(-1)

#UV FUNCTIONS###################################################From MPhys Paper
def E_ion(z):#Hz/erg
    return 10**(25.2)

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

def f_esc_UV(z): #escape fraction of UV
        return 0.2

def n_ion_dot_UV(z): #s⁻¹ cm⁻³ 
    return f_esc_UV(z) * E_ion(z) * P_uv(z) / (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3

def Q_Hii_dot_UV(z,Q_Hii): #s⁻¹	def Q_Hii_dot(z,Q_Hii): #s⁻¹
    return (((n_ion_dot_UV(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16) # conversion from Gyr^-1 to s^-1	    return (((n_ion_dot(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16)  

#MIDDLE STEP FUNCTIONS#################################################First Improvements
def E_ion_LBG():#Hz/erg
    return 10**(25.3)

def n_ion_dot_LBG(z): #s⁻¹ cm⁻³ 
    return f_esc_LBG() * E_ion_LBG() * P_uv(z) / (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3

def f_esc_LBG():
    return 0.03

def Q_Hii_dot_LBG(z,Q_Hii): #s⁻¹	def Q_Hii_dot(z,Q_Hii): #s⁻¹
    return (((n_ion_dot_LBG(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16) # conversion from Gyr^-1 to s^-1	    return (((n_ion_dot(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16)  

#LYA FUNCTIONS#################################################First Improvements

def EW(z): # luminosity density of lyman alpha
    x = pylab.array([2.5,2.8,2.9,3.1, 3.3, 3.7, 4.1, 4.5, 4.8, 5.0, 5.3])
    y = pylab.array([117.134349876, 122.113769531, 172.679885864, 143.5936203, 149.371124268, 96.3648490905, 189.002449036, 181.127731323, 95.0448436864, 125.935153962, 143.343048096])
    p2 = pylab.polyfit(x, y, 1.0)
    p = pylab.poly1d(p2)
    return p(z) 

def E_ion_LAE():
    return 10**25.6		

def f_esc_Lya(z): # esc fraction from Sobral 	
    #return 0.13
    return 0.0048*EW(z)	

def rhoUVfactor(z):
    """
    Function from p. 21 on David's following paper:
    Slicing COSMOS with SC4K: the evolution of typical Lyα emitters and the Lyα escape fraction 
    from z ∼ 2 to z ∼ 6
    """
    x1 = pylab.array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4,17])
    y = pylab.array([10.**(26.52), 10.**(26.30), 10.**(26.10), 10.**(25.98), 10.**(25.67), 10.**(25.32),10.**(23.12) ])
    factors = pylab.array([0.40, 0.30, 0.5, 0.8, 1.0, 1.0,1.0 ]) # Fraction of UV luminosity density recovered by Lyman-alpha emitters
    y = y * factors
    y1 = pylab.array([math.log10(i) for i in y])

    p1 = pylab.polyfit(x1,y1,3.)
    p = pylab.poly1d(p1)
    return 10**p(z)

def n_ion_dot_Lya(z):
    return f_esc_Lya(z)*E_ion_LAE()*rhoUVfactor(z) / (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3

def Q_Hii_dot_Lya(z,Q_Hii): #s⁻¹	
    return (((n_ion_dot_Lya(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16) # conversion from Gyr^-1 to s^-1
