from pylab import *
from scipy.integrate import odeint
from math import *
from scipy.interpolate import interp1d

def AGE(z):
  verbose=0
  length=2
# if no values, assume Benchmark Model, input is z
  if length == 2:
    H0 = 70                         # Hubble constant
    WM = 0.3                        # Omega(matter)
    WV = 1.0 - WM - 0.4165/(H0*H0)  # Omega(vacuum) or lambda
# initialize constants
  WR = 0.        # Omega(radiation)
  WK = 0.        # Omega curvaturve = 1-Omega(total)
  c = 299792.458 # velocity of light in km/sec
  Tyr = 977.8    # coefficent for converting 1/H into Gyr
  DTT = 0.5      # time from z to now in units of 1/H0
  DTT_Gyr = 0.0  # value of DTT in Gyr
  age = 0.5      # age of Universe in units of 1/H0
  age_Gyr = 0.0  # value of age in Gyr
  zage = 0.1     # age of Universe at redshift z in units of 1/H0
  zage_Gyr = 0.0 # value of zage in Gyr
  DCMR = 0.0     # comoving radial distance in units of c/H0
  DCMR_Mpc = 0.0 
  DCMR_Gyr = 0.0
  DA = 0.0       # angular size distance
  DA_Mpc = 0.0
  DA_Gyr = 0.0
  kpc_DA = 0.0
  DL = 0.0       # luminosity distance
  DL_Mpc = 0.0
  DL_Gyr = 0.0   # DL in units of billions of light years
  V_Gpc = 0.0
  a = 1.0        # 1/(1+z), the scale factor of the Universe
  az = 0.5       # 1/(1+z(object))
  h = H0/100.
  WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
  WK = 1-WM-WR-WV
  az = 1.0/(1+1.0*z)
  age = 0.
  n=1000         # number of points in integrals
  for i in range(n):
    a = az*(i+0.5)/n
    adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
    age = age + 1./adot

  zage = az*age/n
  zage_Gyr = (Tyr/H0)*zage
  DTT = 0.0
  DCMR = 0.0
# do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
  for i in range(n):
    a = az+(1-az)*(i+0.5)/n
    adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
    DTT = DTT + 1./adot
    DCMR = DCMR + 1./(a*adot)
  DTT = (1.-az)*DTT/n
  DCMR = (1.-az)*DCMR/n
  age = DTT+zage
  age_Gyr = age*(Tyr/H0)
  DTT_Gyr = (Tyr/H0)*DTT
  DCMR_Gyr = (Tyr/H0)*DCMR
  DCMR_Mpc = (c/H0)*DCMR

# tangential comoving distance
  ratio = 1.00
  x = sqrt(abs(WK))*DCMR
  if x > 0.1:
    if WK > 0:
      ratio =  0.5*(exp(x)-exp(-x))/x 
    else:
      ratio = sin(x)/x
  else:
    y = x*x
    if WK < 0: y = -y
    ratio = 1. + y/6. + y*y/120.
  DCMT = ratio*DCMR
  DA = az*DCMT
  DA_Mpc = (c/H0)*DA
  kpc_DA = DA_Mpc/206.264806
  DA_Gyr = (Tyr/H0)*DA
  DL = DA/(az*az)
  DL_Mpc = (c/H0)*DL
  DL_Gyr = (Tyr/H0)*DL

# comoving volume computation
  ratio = 1.00
  x = sqrt(abs(WK))*DCMR
  if x > 0.1:
    if WK > 0:
      ratio = (0.125*(exp(2.*x)-exp(-2.*x))-x/2.)/(x*x*x/3.)
    else:
      ratio = (x/2. - sin(2.*x)/4.)/(x*x*x/3.)
  else:
    y = x*x
    if WK < 0: y = -y
    ratio = 1. + y/5. + (2./105.)*y*y
  VCM = ratio*DCMR*DCMR*DCMR/3.
  V_Gpc = 4.*pi*((0.001*c/H0)**3)*VCM
  return zage_Gyr

#zl = arange(0.0,15.,0.01) #redshift
t = arange(0.1,3.0,0.1)*pow(10,12)


A=[]
B=[]
for i in arange(0.0,15.,0.1):
    A.append(i)
    B.append(AGE(i))
A=array(A) # list of all redshifts
B=array(B) # list of corresponding times since BB

t_to_z = interp1d(B,A)

xi_ion = 10.**(25.5) #production_efficiency [Hz erg**(-1.)]
f_esc0 = 2.3
alpha = 1.17
x1 = array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4])
y1 = array([10.**(26.52), 10.**(26.30), 10.**(26.10), 10.**(25.98), 10.**(25.67), 10.**(24.62)])
p1 = polyfit(x1,y1,2.)
p = poly1d(p1)
Omega_b = 0.04 #baryon_density
h = 0.7
X_p = 0.75 #hydrogen_mass_fraction
T = 2.*10.**(4.) #IGM_temperature [K]
C = 3 #clumping_factor
 
def dQ_dz(Q,t): #HII_fraction
    z= t_to_z(t)
    return (((f_esc0*((1.+z)/3.)**alpha)/100.)*xi_ion*(p(z)))/((1.67*10.**(-7.)*(Omega_b*h**2/0.02)*(X_p/0.75))/(3.40368*10.**(-74.)))-Q/((C*(3.40368*10.**(-74.)*2.6*10.**(-13.)*(T/10.**(4.))**(-0.76))*(1.+(1.-X_p)/(4.*X_p))*((1.67*10.**(-7.)*(Omega_b*h**2/0.02)*(X_p/0.75))/(3.40368*10.**(-74.)))*(1.+z)**(3.))**(-1.)) #[s**(-1.)]

#def n_ion(f_esc,xi_ion,rho_UV): #production_rate
#    return (((f_esc0*((1.+z)/3.)**alpha)/100)*xi_ion*(p(z))) #[s**(-1.) Mpc**(-3.)]

#def f_esc(f_esc0,alpha,z): #escape_fraction
#    return (f_esc0*((1.+z)/3.)**alpha) #[%]

#def rho_UV(p,z): #UV_luminosity_density 
#    return p(z) #[erg s**(-1.) Hz**(-1.) Mpc**(-3.)]

#def n_H(Omega_b,h,X_p): #hydrogen_number_density
#    return ((1.67*10.**(-7.)*(Omega_b*h**2/0.02)*(X_p/0.75))/(3.40368*10.**(-74.))) #[Mpc**(-3.)]

#def t_rec(C,alpha_B,Y_p,X_p,n_H,z): #hydrogen_recombination_time
#    return ((C*alpha_B*(1.+Y_p/(4.*X_p))*n_H*(1.+z)**(3.))**(-1.)) #[s]

#def t_rec(C,alpha_B,Y_p,X_p,n_H,z): #hydrogen_recombination_time
#    return ((C*alpha_B*(1.+Y_p/(4.*X_p))*n_H*(1.+z)**(3.))**(-1.)) #[s]

#def alpha_B(T): #case_B_recombination_coefficient
#    return (3.40368*10.**(-74.)*2.6*10.**(-13.)*(T/10.**(4.))**(-0.76)) #[Mpc**(3.) s**(-1.)]

#def Y_p(X_p): #helium_mass_fraction
#    return (1.-X_p)




Q = odeint(dQ_dz,1.0,t)
plot(t,Q)
xlabel('Redshift')
ylabel('Q_HII')
title('HII_fraction')
show()
