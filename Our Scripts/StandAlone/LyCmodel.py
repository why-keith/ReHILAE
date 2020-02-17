#Based on pg19 of MPhys Thesis
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
c_ha=1.36E-12 #Recombination coefficient (erg)
EW_avg=140.321828866 #Average equivalent width (angstroms)

#PARAMETER DEFAULTS#############################################
alpha=1.17
T=20000 #K Temperature
C=3 #clumping factor

#FUNCTIONS######################################################
def n_H(): #cm^-3  hydrogen number density
    return 1.67 * (10**(-7)) * (W_b_h_sqr / 0.02) * (X_p / 0.75)

def alpha_beta(): #cm³ s⁻¹ - recombination coefficient
    return 2.6*(10**(-13)) * ((T/(10**4))**(-0.76))

def t_rec(z): #s recombination time
    return (alpha_beta() * n_H() * C * (1 + Y_p/(4*X_p)) * (1 + z)**(3))**(-1)

def P_L_Lya(z):
    """"
    Data from Table_C2_Sobral18_SSC4K_compilation.fits
    """
    x = pylab.array([2.2,2.5,2.8,3.0,3.2, 3.3, 3.7, 4.1, 4.6, 4.8, 5.1, 5.3, 5.8 ])
    y = pylab.array([0.52, 0.74, 0.77, 0.88, 0.84, 0.85, 1.01, 0.87, 1.19, 1.12, 1.27, 1.08, 1.10])  # data from SC4K Sobral 
    p2 = pylab.polyfit(x, y, 2.0)
    p = pylab.poly1d(p2)

    if p(z) < 0:
        return 0
    else:
        return p(z) * 10**40 

def EW(z): # luminosity density of lyman alpha
    """
    Data from "Table_C3_Calhau19_Stacking_LAEs_X_rays_v1.fits"
    """
    x = pylab.array([2.5, 2.8, 2.9, 3.1, 3.3, 3.7, 4.1, 4.5, 4.8, 5.0, 5.3])
    y = pylab.array([117.134349876, 122.113769531, 172.679885864, 143.5936203, 149.371124268, 96.3648490905, 189.002449036, 181.127731323, 95.0448436864, 125.935153962, 143.343048096])
    p2 = pylab.polyfit(x, y, 1.0)
    p = pylab.poly1d(p2)
    return p(z) 

def Q_ion_LyC(z): # replaces P_uv and E_ion	
    return P_L_Lya(z) / (c_ha*(1 - f_esc_LyC(EW(z)))*(0.042 * EW(z)))	

def f_esc_LyC(z):  #escape fraction of lyman continuum from Lya
    """
    Data from https://www.aanda.org/articles/aa/pdf/2017/01/aa29264-16.pdf
    """
    x = pylab.array([79, 129, 83, 98, 75, 29, 15, 4])
    y = pylab.array([0.132, 0.074, 0.072, 0.058, 0.056, 0.045, 0.032, 0.01])
    p1 = pylab.polyfit(x,y,1.0)
    p = pylab.poly1d(p1)
    return p(EW(z))


def n_ion_dot_LyC(z): # replaces n_ion_dot using Q_ion_LyC	
    return Q_ion_LyC(z) * f_esc_LyC(z) / (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3
    
def Q_Hii_dot(z,Q_Hii): #s⁻¹	def Q_Hii_dot(z,Q_Hii): #s⁻¹
    return (((n_ion_dot_LyC(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16) # conversion from Gyr^-1 to s^-1	    return (((n_ion_dot(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16)