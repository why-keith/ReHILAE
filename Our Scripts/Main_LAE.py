import math
import numpy as np
import pylab
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import curve_fit
from scipy.integrate import odeint

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


#log_list=[0,0,0,0]
init_conditions={"P1":np.array([]), "P2":np.array([]), "P3":np.array([]), "f1":np.array([]), "f2":np.array([])}
#log=open(r"log.txt","a")

P1 = -0.05
P2 = 0.44
P3 = 38.99

f1 = 0.00064
f2 = 0.00941

#PARAMETER DEFAULTS#############################################
alpha=1.17
T=20000 #K Temperature
C=3 #clumping factor
startT=0.124 #Start time
finishT=14 #Finish time
intervalNumber = 10000
TStep=(finishT - startT)/(intervalNumber) #Size of steps in time
EW = 148.9705


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
    Z=red(T, alt)
    return Z,T

z,t = redshift(startT, finishT, TStep)


def alpha_beta(): #cm³ s⁻¹ - recombination coefficient
    return 2.6*(10**(-13)) * ((T/(10**4))**(-0.76))

def n_H(): #cm^-3  hydrogen number density
    return 1.67 * (10**(-7)) * (W_b_h_sqr / 0.02) * (X_p / 0.75)

def t_rec(z): #s recombination time
    return (alpha_beta() * n_H() * C * (1 + Y_p/(4*X_p)) * (1 + z)**(3))**(-1)

def P_L_Lya(x, P1, P2, P3):

    return 10**(P1*x**2 + P2*x +P3)


def f_esc_LyC(x, f1, f2):  #////////////////////////////////////////////////////////////////////////////////// LOGGED
   # log.write("f_esc_LyC\n")
    #print("f_esc_LyC")
   # log_list[0]+=1
    
    f_esc=f1*x + f2
    return f_esc

def Q_ion_LyC(z, P1, P2, P3, f1, f2): # replaces P_uv and E_ion	

    return P_L_Lya(z, P1, P2, P3) / ((c_ha*(1 - f_esc_LyC(EW, f1, f2)))*(0.042 * EW))


def n_ion_dot_LyC(z, P1, P2, P3,  f1, f2): # replaces n_ion_dot using Q_ion_LyC	 ////////////////////////// LOGGED
#    log.write("n_ion_dot    z = "+str(z)+"\n")
   # print("n_ion_dot")
    #log_list[1]+=1
    if Q_ion_LyC(z, P1, P2, P3,  f1, f2) * f_esc_LyC(EW, f1, f2) / (2.938e+73) <= 0 :
        n_ion_dot= 0 
    else:
        n_ion_dot= Q_ion_LyC(z, P1, P2, P3, f1, f2) * f_esc_LyC(EW, f1, f2) / (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3

    return n_ion_dot

def Q_Hii_dot(z, Q, P1, P2, P3, f1, f2): #s⁻¹	def Q_Hii_dot(z,Q_Hii): #s⁻¹    /////////////////////////////////////////
#    log.write("Q_Hii_dot    z = "+str(z)+"\n")
   # print("Q_Hii_dot")
    #log_list[2]+=1
    Q_dot=(((n_ion_dot_LyC(z, P1, P2, P3, f1, f2)/n_H()) - (Q/t_rec(z)))*3.1536e+16)    
    
    return Q_dot # conversion from Gyr^-1 to s^-1	    return (((n_ion_dot(z)/n_H()) - (Q_Hii/t_rec(z)))*3.1536e+16)  

def dQ_dt(Q, t, P1, P2, P3, f1, f2):    
    dQ_dt = Q_Hii_dot(red(t), Q, P1, P2, P3, f1, f2)
    return dQ_dt

#GENERATE Q ARRAY

def main(arguements):
    Q = odeint(dQ_dt, Q_Hii_zero, t, args=(arguements))
 #   print("arguements = {}".format(arguements))
 #   print(type(arguements))
    Q[Q>1.0] = 1.0 # 100% HII
    Q[Q<0.0] = 0.0 # 100% HI
   # print("Q")
  #  log.write("Q+\n")
    #log_list[3]+=1
    
    for i,j in zip(init_conditions,range(len(arguements))):        
    #    print("appending {} to {}".format(j,i))
        init_conditions[i]=np.append(init_conditions[i],arguements[j])
     #   print(init_conditions[i])

    
    
   # print("----------------------------")
    return Q    #///////////////////////////////////////////////////////////////////////

