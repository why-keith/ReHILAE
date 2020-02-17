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
startT=0.124 #Start Z
finishT=14 #Finish z
intervalNumber = 10000
TStep=(finishT - startT)/(intervalNumber) #Size of steps in Z



#FUNCTIONS######################################################
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

def f_esc_UV(z): #escape fraction of UV
        return (f_esc_zero*((1+z)/3)**alpha)/100 # UV escape fraction

def n_ion_dot_UV(z): #s⁻¹ cm⁻³ 
    return f_esc_UV(z) * E_ion(z) * P_uv(z) / (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3

def Q_Hii_dot_UV(z,Q_Hii): #s⁻¹	def Q_Hii_dot(z,Q_Hii): #s⁻¹
    return (((n_ion_dot_UV(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16) # conversion from Gyr^-1 to s^-1	    return (((n_ion_dot(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16)  


    
#LYA FUNCTIONS#################################################First Improvements
def P_L_Lya(z): # luminosity density of lyman alpha
    # data from SC4K Sobral
    x = pylab.array([2.2,2.5,2.8,3.0,3.2, 3.3, 3.7, 4.1, 4.6, 4.8, 5.1, 5.3, 5.8 ])
    y = pylab.array([0.52, 0.74, 0.77, 0.88, 0.84, 0.85, 1.01, 0.87, 1.19, 1.12, 1.27, 1.08, 1.10])  # data from SC4K Sobral 
    p2 = pylab.polyfit(x, y, 2.0)
    p = pylab.poly1d(p2)
    
    peak=np.max(y)
    peak_position = np.where(y==peak)
    cutoff=x[peak_position[0][0]]
    
    if z > cutoff:
        return p(cutoff) * 10**40
    if p(z) < 0:
        return 0
    else:
        return p(z) * 10**40 

#def f_esc_Lya(z): # esc fraction from Sobral 	
#    return 0.0048*EW(z)	

def EW(z): # luminosity density of lyman alpha
    x = pylab.array([2.5,2.8,2.9,3.1, 3.3, 3.7, 4.1, 4.5, 4.8, 5.0, 5.3])
    y = pylab.array([117.134349876, 122.113769531, 172.679885864, 143.5936203, 149.371124268, 96.3648490905, 189.002449036, 181.127731323, 95.0448436864, 125.935153962, 143.343048096])
#    filename = 'Table_C3_Calhau19_Stacking_LAEs_X_rays_v1.fits'
#    hdu_list = fits.open(filename) 
#    evt_data = Table(hdu_list[1].data) 
#    stack = evt_data.field(0)[4:15]
#    redshift = []
#    for i in stack: # cleaning data e.g. 'z=2.5-NO_AGN' => 2.5
#        j = str(i).split('=')
#        j = j[1]
#        j = j.split('-')
#        j = j[0]
#        redshift.append(float(j))

#    EquiWidth = evt_data.field(8)[4:15]
#    EW = [float(i) for i in EquiWidth]
#    x, y = pylab.array(redshift), pylab.array(EW)
    p2 = pylab.polyfit(x, y, 1.0)
    p = pylab.poly1d(p2)
    return p(z) 



#LYC FUNCTIONS#################################################Next Improvements
def f_esc_Lya(z): # esc fravction from Sobral 	
    return 0.0048*EW(z)	

def Q_ion_LyC(z): # replaces P_uv and E_ion	
    return P_L_Lya(z) / (c_ha*(1 - f_esc_LyC(EW(z)))*(0.042 * EW(z)))	

def f_esc_LyC(z):  #escape fraction of lyman continuum from Lya
    x = pylab.array([79, 129, 83, 98, 75, 29, 15, 4])
    y = pylab.array([0.132, 0.074, 0.072, 0.058, 0.056, 0.045, 0.032, 0.01])
    p1 = pylab.polyfit(x,y,1.0)
    p = pylab.poly1d(p1)
    return p(EW(z))

def n_ion_dot_LyC(z): # replaces n_ion_dot using Q_ion_LyC	
    if Q_ion_LyC(z) * f_esc_LyC(z) / (2.938e+73) < 0 :
        return 0 
    else:
        return Q_ion_LyC(z) * f_esc_LyC(z) / (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3
    
def Q_Hii_dot(z,Q_Hii): #s⁻¹	def Q_Hii_dot(z,Q_Hii): #s⁻¹
    return (((n_ion_dot_LyC(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16) # conversion from Gyr^-1 to s^-1	    return (((n_ion_dot(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16)  



