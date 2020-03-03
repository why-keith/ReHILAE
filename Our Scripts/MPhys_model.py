"""
Defines, calculates and allows for the changing of variables used in the toy model 
"""
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

"""	
def set_variables(_alpha=alpha,_T=T,_C=C,_startT=startT,_finishT=finishT,_intervalNumber=intervalNumber): #allows changing of parameters from outside model.py	
    global alpha	
    global T	
    global C	
    global startT	
    global finishT	
    global intervalNumber	
    global TStep	

    alpha=_alpha	
    T=_T	
    C=_C	
    startT=_startT	
    finishT=_finishT	
    intervalNumber=_intervalNumber	
    TStep = (finishT - startT)/(intervalNumber)	


    return "α={} \nT={}\nC={}".format(alpha,T,C)
"""  
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
    return 10**25.6 #10**(24.4 + math.log10(1 + z))

def squared(x, a1, a2, a3):
    return a1*x**2 + a2*x +a3

def P_uv(z):   #Hz^-1 s^-1 Mpc^-3  UV luminosity density 
    """
    From Bouwens 2015 Table 7

    available at: https://arxiv.org/pdf/1403.4295.pdf

    returns log10(P_uv)
    """
    z = math.log10(1+z)
    return -5.288*z + 30.39894


def f_esc_UV(z): #escape fraction of UV
        return (f_esc_zero*((1+z)/3)**alpha)/100 # UV escape fraction

def n_ion_dot_UV(z): #s⁻¹ cm⁻³ 
    return f_esc_UV(z) * E_ion(z) * (10**P_uv(z)) / (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3

def Q_Hii_dot_UV(z,Q_Hii): #s⁻¹	def Q_Hii_dot(z,Q_Hii): #s⁻¹
    return (((n_ion_dot_UV(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16) # conversion from Gyr^-1 to s^-1	    return (((n_ion_dot(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16)  


    
#LYA FUNCTIONS#################################################First Improvements
def P_L_Lya(z): # luminosity density of lyman alpha
    """
    SC4K data shown in Sobral et al. 2018 table C3

    available at: https://arxiv.org/pdf/1712.04451.pdf

    returns P_L_Lya
    """

    z = math.log10(1+z)
    return -2.97564*z + 41.96476



def linear(x, a1, a2):
    return a1*x + a2

def EW(z):
    """
    SC4K data shown in Calhau et al. 2019 Table C3

    available at: https://lancaster.app.box.com/s/t75t3v713yuibkvjk3ioqdesunxpv1fb/file/618125008170
    """
    x = pylab.array([2.5,2.8,2.9,3.1, 3.3, 3.7, 4.1, 4.5, 4.8, 5.0, 5.3, 5.8])
    y = np.array([117.134349876, 122.113769531, 172.679885864, 143.5936203, 149.371124268, 96.3648490905, 189.002449036, 181.127731323, 95.0448436864, 125.935153962, 143.343048096, 251.93561])
    y_err_up =pylab.array([186.35171, 160.96419, 321.31418, 382.54525, 258.18398, 70.23115, 763.57573, 223.72103, 100.98794, 195.63354, 159.8391, 620.82504])
    y_err_down =pylab.array([53.82247, 53.55149, 90.78461, 72.25881, 74.65443, 39.08451, 96.45726, 85.27828, 42.09321, 52.08989, 71.00167, 186.34499])
    error = [(y_err_up[i]-y_err_down[i])/2 for i in range(len(y_err_down))]
    p2 = np.polyfit(x, y, 1.0)
    p = pylab.poly1d(p2)
    Params = p.coefficients
    pfit, pcov = curve_fit(linear, x, y, p0=Params, sigma=error)
    a1,a2= pfit[0], pfit[1]
    aver = np.mean(y)
    print(aver)
    return linear(z, a1, a2)



#LYC FUNCTIONS#################################################Next Improvements
def f_esc_Lya(z): # esc fravction from Sobral 	
    return 0.0048*EW(z)	

def Q_ion_LyC(z): # replaces P_uv and E_ion	
    return 10**P_L_Lya(z) / (c_ha*(1 - f_esc_LyC(EW(z)))*(0.042 * EW(z)))

def f_esc_LyC(z):  #escape fraction of lyman continuum from Lya
    """
    SDSS data shown in Verhamme et al. 2016 Table 1

    availabke at: https://www.aanda.org/articles/aa/pdf/2017/01/aa29264-16.pdf
    """
    x = pylab.array([79, 129, 83, 98, 75, 29, 15, 4]) # Equivalent Width
    y = pylab.array([0.132, 0.074, 0.072, 0.058, 0.056, 0.045, 0.032, 0.01]) # f_esc
    y_sigma = [0.2*y[i] for i in range(len(y))]
    p1 = np.polyfit(x,y,1.0)
    p = pylab.poly1d(p1)
    Params = p.coefficients
    pfit, pcov = curve_fit(linear, x, y, p0=Params, sigma=y_sigma)
    a1,a2= pfit[0], pfit[1]
    return linear(z, a1, a2)

def n_ion_dot_LyC(z): # replaces n_ion_dot using Q_ion_LyC	
    if Q_ion_LyC(z) * f_esc_LyC(EW(z)) / (2.938e+73) <= 0 :
        return 0 
    else:
        return Q_ion_LyC(z) * f_esc_LyC(EW(z)) / (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3
    
def Q_Hii_dot(z,Q_Hii): #s⁻¹	def Q_Hii_dot(z,Q_Hii): #s⁻¹
    return (((n_ion_dot_LyC(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16) # conversion from Gyr^-1 to s^-1	    return (((n_ion_dot(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16)  
