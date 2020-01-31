from pylab import *
from scipy.integrate import odeint

f_esc0 = 2.3
alpha = 1.17
xi_ionc = 10.**(25.5) # Constant xi_ion
T = 2.*10.**(4.) # IGM Temperature [K]
Omega_b = 0.04 # Baryon Density
h = 0.7
X_p = 0.75 # Hydrogen Mass Fraction
C = 3 # Clumping Factor
x = array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4])
y = array([10.**(26.52), 10.**(26.30), 10.**(26.10), 10.**(25.98), 10.**(25.67), 10.**(24.62)])
p1 = polyfit(x,y,2.)
p = poly1d(p1)
t = arange(0.01,14.,0.1)

#def z(t):
#    return ((28./t)-1.)**(1./2.)-1.

#def f_esc(f_esc0,alpha,z): # Escape Fraction [%]
#    return f_esc0*((1.+z)/3.)**alpha

#def rho_UV(p,z): # UV Luminosity Density [erg s^-1 Hz^-1 Mpc^-3]
#    return p(z)

#def n_ion(f_esc,xi_ionc,rho_UV): # Production Rate (Constant Production Efficiency) [s^-1 Mpc^-3]
#    return (f_esc/100)*xi_ionc*rho_UV

#def alpha_B(T): # Case B Recombination Coefficient [Mpc^3 s^-1]
#    return (3.40368*10.**(-74.))*(2.6*10.**(-13.))*((T/(10.**(4.)))**(-0.76))

#def n_H(Omega_b,h,X_p): # Hydrogen Number Density [Mpc^-3]
#    return ((1.67*10.**(-7.))*((Omega_b*(h**(2.)))/0.02)*(X_p/0.75))/(3.40368*10.**(-74.))

#def t_rec(C,alpha_B,X_p,n_H,z): # Hydrogen Recombination Time [s]
#    return (C*alpha_B*(1.+((1.-X_p)/(4.*X_p)))*n_H*((1.+z)**(3.)))**(-1.)

#def dQ_dt(Q,t): #HII Fraction
#    return (n_ion/n_H)-(Q/t_rec)

def dQ_dt(Q,t): # HII Fraction
    return ((((f_esc0*((1.+(((28./t)-1.)**(1./2.)-1.))/3.)**alpha)/100)*xi_ionc*(p(((28./t)-1.)**(1./2.)-1.)))/(((1.67*10.**(-7.))*((Omega_b*(h**(2.)))/0.02)*(X_p/0.75))/(3.40368*10.**(-74.))))-(Q/((C*((3.40368*10.**(-74.))*(2.6*10.**(-13.))*((T/(10.**(4.)))**(-0.76)))*(1.+((1.-X_p)/(4.*X_p)))*(((1.67*10.**(-7.))*((Omega_b*(h**(2.)))/0.02)*(X_p/0.75))/(3.40368*10.**(-74.)))*((1.+(((28./t)-1.)**(1./2.)-1.))**(3.)))**(-1.)))

Q = odeint(dQ_dt,0.,t) # the 2nd argument encodes Q(0) = 0
plot(t,Q)
xlabel('Cosmic Time (t) [Gyr]')
ylabel('Q_HII')
title('HII Fraction')
show()
