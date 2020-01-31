#
# Solving Re-ionisation
# By DS 2017
#
from pylab import *
from scipy.integrate import odeint
import matplotlib
matplotlib.rc('text', usetex=True)
from math import *
from scipy.interpolate import interp1d
import pyfits as py
from scipy.integrate import quadpack
from scipy.integrate import quadrature
from scipy import special
import os

##################################
def dQ_dt (Q,t): #HII_fraction

    z = t_to_z(t)
    trec = 1./(C*alpha_B*n_e*(1.+z)**(3.)) # s 
    nion = f_esc0*xi_ion*p(z)/ (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3
    return ((nion)/n_H-Q/trec)*3.1536e+16 # Conversion to time in seconds
    
##################################
def tau (z):   # NEED final check

    return sigma_T*n_e*(1.+z)*(1.+z)*Q(z)*c/(H)/z   # Units ok - in

    # sigma_T  = cm^2 
    # c km/s
    # H km/s/Mpc => now km/s/cm
    
    # Q(z) - unitless
    #n_e per cm^3
    
    # cm^2 * cm  / cm^3 = unitless

def compute_tau_z(Q,zmin=0.6,zmax=14.,dz=0.3):
    z_integrate = arange(zmin,zmax,dz)
    taus_plot = []
    for i in z_integrate:
        taus_plot.append(quadpack.quad(tau,0.3,i)[0])
        
    return z_integrate, taus_plot   
    
###########################
# Main properties/constants
#####################################################################################
Omega_b = 0.04 #baryon_density
h = 0.7 # hubble constant / 100
X_p = 0.75 #hydrogen_mass_fraction
T = 20000. #IGM_temperature [K]
C = 3.  #clumping_factor
Y_p = (1.-X_p) # Helium mass
alpha_B = 2.6*pow(10,-13.)*(T/10.**(4.))**(-0.76)  #case_B_recombination_coefficient
n_H = 1.67*pow(10,-7.) * ((Omega_b * h* h)/0.02) * (X_p/0.75) # per cm^3

c = 299792.458

#### For getting Tau (z)

sigma_T = 6.65*pow(10,-25) # Thomson scattering cross section   cm^2  - from Mayer et al. 
n_e = (1.+Y_p/(4.*X_p))*n_H   # comoving number density of electrons (assuming singly ionized He
H = h*100. *(3.24077929*pow(10,-25.)) # Hubble parameter == hubble constant conversion per cm from per Mpc

#######################################################################################

##################################################################
# GET redshift - time conversion
######
# F737 Cosmology
F737 = 'Redshift_time_F737_cosmology.fits'
# PLanck cosmology
FPlanck = 'Redshift_time_FPlanck_cosmology.fits'
#LOAD Redshift-time dependence
A=py.open(F737) 
NOW=A[1].data
#LOAD GENERAL PROPERTIES - what will be plotted
Redshift = NOW.field('Redshift') 
Time =  NOW.field('Time_Gyr') 
# time to redshift relation
t_to_z = interp1d(Time,Redshift,kind='linear',bounds_error=None)

########################################
# UV luminosity densities
B=py.open('RHO_UV.fits') 
NOWB=B[1].data
Red_z = NOWB.field('Redshift') 
rho_UV =  10.**(NOWB.field('p_UV_13')) # 'p_UV_12  p_UV_13', 'p_UV_17' - other integration limits
p = interp1d(Red_z,rho_UV,kind='linear',bounds_error=None) # p is rho UV(z)

############################################################
# Integration limits and conversion z to time
#
ts = np.linspace(0.25,13.4,20000000)  # Produce time grid
zs= t_to_z(ts)
t0 = 0.0
###############################################################

#### START PLOT

# DEFINE the sub-plots to make them side by side and share Y axis
f, (ax1,ax2) = plt.subplots(2,1, sharex=True, figsize=(9,9))
plt.subplots_adjust(wspace = -0.04)  # adjust spacing between sub-figures
# setting up the font style throughout 
matplotlib.rcParams.update({'font.size': 18, 'font.weight': 'bold'}) # change for diff style

###### QHII PLOT NOW ######
########################################
# Plot end of re-ionisation
za = [6.0,6.0000001]
zp = [-0.5,5.0]
Reionisation = plot(za,zp,'k--')

## PLOT REIONISATION EPOCH
za = [8.8,8.80001]
zp = [-0.5,5.0]
Reionisation = plot(za,zp,linestyle='-',color ='0.40', linewidth=8, zorder=-8)
Reionisation = plot(array(za)+0.05,zp,linestyle='-',color ='0.90', linewidth=100, zorder=-10)
Reionisation = plot(array(za)+0.05,zp,linestyle='-',color ='0.97', linewidth=200, zorder=-20)

Reionisation = ax1.plot(za,zp,linestyle='-',color ='0.60', linewidth=8, zorder=-18)
#Reionisation = ax1.plot(array(za)+0.05,zp,linestyle='-',color ='0.94', linewidth=100, zorder=-10)
#Reionisation = ax1.plot(array(za)+0.05,zp,linestyle='-',color ='0.99', linewidth=200, zorder=-20)

#################################################################
# ALL - Robertson+ like
xi_ion = 10.**(25.2) #production_efficiency [Hz erg**(-1.)]
f_esc0 = 0.20

Q = odeint(dQ_dt,t0,ts)  # SOLVE NOW
Q_LAE=Q
Q[Q>1.0] = 1.0 # boundary condition up
Q[Q<0.0] = 0.0 # boundary condition down

# Make it an actual function - rder values
Q_not = []
z_not = []
zs = array(zs)
for i in range(len(Q)):
    Q_not.append(Q[len(Q)-i-1][0])
    z_not.append(zs[len(zs)-i-1])    

# now interpolate
Q = interp1d(array(z_not),array(Q_not),kind='linear') # p is rho UV(z)
TOTAL = plot(z_not,Q(z_not),'g-',linewidth=3,zorder=-1)

##################################################################################
### NOW TAU PLOT

# SOLVE TAU:
z_plot, taus_plot = compute_tau_z(Q)
Prediction = ax1.plot(z_plot, taus_plot,'g-',linewidth=3,zorder=-1)

#################################################################
# LAEs
xi_ion = 10.**(25.6) #production_efficiency [Hz erg**(-1.)]
f_esc0 = 0.10

Q = odeint(dQ_dt,t0,ts)  # SOLVE NOW
Q_LAE=Q
Q[Q>1.0] = 1.0 # boundary condition up
Q[Q<0.0] = 0.0 # boundary condition down

# Make it an actual function - rder values
Q_not = []
z_not = []
zs = array(zs)
for i in range(len(Q)):
    Q_not.append(Q[len(Q)-i-1][0])
    z_not.append(zs[len(zs)-i-1])    

# now interpolate
Q = interp1d(array(z_not),array(Q_not),kind='linear') # p is rho UV(z)
LAEs = plot(z_not,Q(z_not),'b-.',linewidth=6,zorder=-1)

##################################################################################
### NOW TAU PLOT

# SOLVE TAU:
z_plot, taus_plot = compute_tau_z(Q)
Prediction = ax1.plot(z_plot, taus_plot,'b-.',linewidth=6,zorder=-1)



#################################################################
# LBGs z~3
xi_ion = 10.**(25.3) #production_efficiency [Hz erg**(-1.)]
f_esc0 = 0.03

Q = odeint(dQ_dt,t0,ts)  # SOLVE NOW
Q_LAE=Q
Q[Q>1.0] = 1.0 # boundary condition up
Q[Q<0.0] = 0.0 # boundary condition down

# Make it an actual function - rder values
Q_not = []
z_not = []
zs = array(zs)
for i in range(len(Q)):
    Q_not.append(Q[len(Q)-i-1][0])
    z_not.append(zs[len(zs)-i-1])    

# now interpolate
Q = interp1d(array(z_not),array(Q_not),kind='linear') # p is rho UV(z)
LBGs_real = plot(z_not,Q(z_not),'r--',linewidth=3,zorder=-1)

##################################################################################
### NOW TAU PLOT

# SOLVE TAU:
z_plot, taus_plot = compute_tau_z(Q)
Prediction = ax1.plot(z_plot, taus_plot,'r--',linewidth=3,zorder=-1)





#### PLOT Planck
########################################
# Plot Planck optical depth
za = [-2.0,19.0000001]
zp = [0.066,0.06600001]
Reionisation = ax1.plot(za,zp,'k--',linewidth=3)
Reionisation = ax1.plot(za,zp,'k-',linestyle='-',color ='0.90', linewidth=50, zorder=-20) # INCLUDES ERROR
#Reionisation = ax1.plot(za,array(zp)+0.012,'k-',linewidth=3)
#Reionisation = ax1.plot(za,array(zp)-0.012,'k-',linewidth=3)


# LEGEND
I1=legend((TOTAL[0],LAEs[0],LBGs_real[0]),(r'R+13,15',r'LAEs $z\sim2-3$',r'LBGs $z\sim2-3$'), shadow = True, loc=3,numpoints=1)
ltext = gca().get_legend().get_texts()
setp(ltext[0], fontsize = 19, color='0.1')
setp(ltext[1], fontsize = 19, color='k')
frame=I1.get_frame()
frame.set_linewidth(0)
frame.set_visible(False)

# Define the axis
axis([0,13.1,-0.05,1.05])
ax1.axis([0,13.1,-0.01,0.11])

# Label axis - quantities and units
ax2.set_xlabel(r'Redshift ($\bf z$)', {'color':'k', 'fontsize': 23})
ax2.set_ylabel(r'Ionised H Fraction (\%)', {'color':'k', 'fontsize': 23})
ax1.set_ylabel(r'$\bf \tau (z)$', {'color'    : 'k', 'fontsize'   :24 })

# SET AXIS - labelling
ax2.set_xticks((0,2,4,6,8,10,12),(0,2,4,6,8,10,12))
ax2.set_yticks((0.0,0.2,0.4,0.6,0.8,1.0),(0,20,40,60,80,100))
ax1.set_yticks((0.0,0.02,0.04,0.06,0.08,0.10),(0.0,0.02,0.04,0.06,0.08,0.10))

savefig('Q_HII_and_Tau_evolution.ps')
show()
