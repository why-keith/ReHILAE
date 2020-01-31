#############################################
#
# Integrating UV luminosity functions
#
#
#############################################

# Import necessary libraries
import os, sys
from numpy import *
import scipy
from scipy.integrate import quadpack
from scipy.integrate import quadrature
from scipy import special

# Schechter function - linear version
#
def Schechter(x):
        p=ini
        return (p[0]/p[1])*(x/p[1])**(p[2])*exp(-x/p[1])*x # 

# From Mason et al. 2015
Redshifts =[0,2,4,5,6,7,8,9,10,12,14,16]
alphas = [-1.68,-1.46,-1.64,-1.75,-1.83,-1.95,-2.10,-2.26,-2.47,-2.74,-3.11,-3.51]
Mstar = [-19.9,-20.3,-21.2,-21.2,-20.9,-21.0,-21.3,-21.2,-21.1,-21.0,-20.9,-20.7]
Phistar = [-2.97,-2.52,-2.93,-3.12,-3.19,-3.48,-4.03,-4.50,-5.12,-5.94,-7.05,-8.25]

Redshifts = [1.65,2.0,2.7]
alphas = [-1.56,-1.72,-1.94]
Mstar = [-19.74,-20.41,-20.71]
Phistar = [log10(2.32)-3.,log10(1.50)-3.,log10(0.55)-3.]

# Flux from M_UV will be at 10 pc, meaning need to convert 100 pc^2 to cm^2  and have 4 * pi factor for flux to luminosity
FACTOR_LUM = (4.*pi*9.521e+38) # 

Min_int0 =  pow(10,-0.4*(-10+48.6)) * FACTOR_LUM   # Lower limit luminosity for integration 1 (M_UV=-12)
Min_int1 =  pow(10,-0.4*(-13+48.6)) * FACTOR_LUM   # Lower limit luminosity for integration 1 (M_UV=-13)
Min_int2 =  pow(10,-0.4*(-17+48.6)) * FACTOR_LUM   # Lower limit luminosity for integration 2 (M_UV=-17)
Max_int = pow(10,-0.4*(-25.5+48.6)) * FACTOR_LUM     # Defining an upper limit for integration just for numeric resaons

os.system('rm z_puv.cat') # delete file to begin with

print >> open('z_puv.cat','a'), '#  Redshift   p_UV_TOTAL   p_UV_10   p_UV_13   p_UV_17'

# Now go through the redshifts:
#
for i in range(len(Redshifts)):
    UVstar = pow(10,-0.4*(Mstar[i]+48.6))  * FACTOR_LUM  # Flux from M_UV will be at 10 pc, meaning need to convert 100 pc^2 to cm^2  and have 4 * pi factor for flux to luminosity
    
    ini = [pow(10.,Phistar[i]), UVstar,alphas[i]]
    
    #print ini
    
    integral1D_m0=quadpack.quad(Schechter,Min_int0,Max_int)
    integral1D_m=quadpack.quad(Schechter,Min_int1,Max_int)
    integral1D_m2=quadpack.quad(Schechter,Min_int2,Max_int)
    integral_analytical = pow(10.,Phistar[i])*UVstar*special.gamma(2+ alphas[i])
    #integral2D_m=quadpack.quad(Schechter,pow(10,-0.4*(-Max_int+48.6)),pow(10,-0.4*(-30+48.6)))
    print
    print Redshifts[i], log10(integral1D_m[0])
    print >> open('z_puv.cat','a'), Redshifts[i], integral_analytical, integral1D_m0[0], integral1D_m[0], integral1D_m2[0]
    #print integral2D_m