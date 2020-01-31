#
# Solving Re-ionisation
# By DS 2017
#
from pylab import *
from scipy.integrate import odeint
import matplotlib
#matplotlib.rc('text', usetex=True)
from math import *
from scipy.interpolate import interp1d

#
### Main properties/constants
#
xi_ion = 10.**(25.6) #production_efficiency [Hz erg**(-1.)]
f_esc0 = 0.2
alpha = 1.17 # this is for f_esc dependence on redshift
x1 = array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4,17])
y1 = array([10.**(26.52), 10.**(26.30), 10.**(26.10), 10.**(25.98), 10.**(25.67), 10.**(24.62),10.**(23.12) ])
p1 = polyfit(x1,y1,3.)
p = poly1d(p1)
Omega_b = 0.04 #baryon_density
h = 0.7
X_p = 0.75 #hydrogen_mass_fraction
T = 20000. #IGM_temperature [K]
C = 3  #clumping_factor
Y_p = (1.-X_p)

#alpha_B_JAKE = (3.40368*10.**(-74.)*2.6*10.**(-13.)*(T/10.**(4.))**(-0.76))  #case_B_recombination_coefficient
alpha_B = 2.6*pow(10,-13.)*(T/10.**(4.))**(-0.76)  #case_B_recombination_coefficient
#n_H_JAKE = ((1.67*10.**(-7.)*(Omega_b*h**2/0.02)*(X_p/0.75))/(3.40368*10.**(-74.)))
n_H = 1.67*pow(10,-7.) * ((Omega_b * h* h)/0.02) * (X_p/0.75) # per cm^3
#n_ion = f_esc0*xi_ion*10.**(21.10)/ (2.938e+73)

# Units! n_H given in cm^-3, but nion usually in Mpc^-3. 1 Mpc^3 = 2.938*pow(10,73)

def dQ_dt (Q,t): #HII_fraction

    z= ((((28./(t))-1.)**(1./2.)-1.))
    
    trec = 1./(C*alpha_B*(1.+Y_p/(4.*X_p))*n_H*(1.+z)**(3.)) # s 
    nion = f_esc0*xi_ion*p(z)/ (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3
    return ((nion)/n_H-Q/trec)*3.1536e+16

#def n_ion(z): #production_rate
#    return (((f_esc0*((1.+z)/3.)**alpha)/100)*xi_ion*(p(z))) #[s**(-1.) Mpc**(-3.)]

#def n_ion(z): #production_rate
#    return f_esc0*xi_ion*10.**(27.10) #[s**(-1.) Mpc**(-3.)]

#def f_esc(z,f_esc0=2.3,alpha=1.17): #escape_fraction
#    return (f_esc0*((1.+z)/3.)**alpha) #[%]

#def rho_UV(p,z): #UV_luminosity_density 
#    return p(z) #[erg s**(-1.) Hz**(-1.) Mpc**(-3.)]

#def n_H(Omega_b,h,X_p): #hydrogen_number_density
#    return ((1.67*10.**(-7.)*(Omega_b*h**2/0.02)*(X_p/0.75))/(3.40368*10.**(-74.))) #[Mpc**(-3.)]

#def t_rec(C,alpha_B,Y_p,X_p,n_H,z): #hydrogen_recombination_time
#    return ((C*alpha_B*(1.+Y_p/(4.*X_p))*n_H*(1.+z)**(3.))**(-1.)) #[s]

#def alpha_B(T): #case_B_recombination_coefficient
#    return (3.40368*10.**(-74.)*2.6*10.**(-13.)*(T/10.**(4.))**(-0.76)) #[Mpc**(3.) s**(-1.)]

#def Y_p(X_p): #helium_mass_fraction
#    return (1.-X_p)

#zs = np.logspace(0.1,16,1000000)
#z0 = 1.0
#Q = odeint(dQ_dt,z0,zs)

ts = np.linspace(0.051,14,10000000)
zs= ((((28./(ts))-1.)**(1./2.)-1.))
t0 = 0.0

#for f_esc0 in arange(0.05,0.3,0.05):
for f_esc0 in arange(0.01,0.2,0.04):
    Q = odeint(dQ_dt,t0,ts)

    Q[Q>1.0] = 1.0
    Q[Q<0.0] =0.0
    plot(zs,Q)


# Define the axis
axis([0,11,-0.05,1.05])

# Label axis
xlabel(r'Redshift', {'color':'k', 'fontsize': 21})
ylabel(r'Ionised Fraction (%)', {'color':'k', 'fontsize': 21})

xticks((0,2,4,6,8,10,12),(0,2,4,6,8,10,12), color='k', size=16)
yticks((0.0,0.2,0.4,0.6,0.8,1.0),(0,20,40,60,80,100), color='k', size=16)

savefig('ionised_fraction.pdf')
show()
