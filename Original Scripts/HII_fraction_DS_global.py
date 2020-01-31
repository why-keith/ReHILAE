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

##################################
def dQ_dt (Q,t): #HII_fraction

    z = t_to_z(t)
    trec = 1./(C*alpha_B*(1.+Y_p/(4.*X_p))*n_H*(1.+z)**(3.)) # s 
    nion = f_esc0*xi_ion*p(z)/ (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3
    return ((nion)/n_H-Q/trec)*3.1536e+16 # Conversion to time in seconds
    
###########################
# Main properties/constants
#####################################################################################
Omega_b = 0.04 #baryon_density
h = 0.7 # hubble constant / 100
X_p = 0.75 #hydrogen_mass_fraction
T = 20000. #IGM_temperature [K]
C = 3  #clumping_factor
Y_p = (1.-X_p) # Helium mass
alpha_B = 2.6*pow(10,-13.)*(T/10.**(4.))**(-0.76)  #case_B_recombination_coefficient
n_H = 1.67*pow(10,-7.) * ((Omega_b * h* h)/0.02) * (X_p/0.75) # per cm^3
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
rho_UV =  10.**(NOWB.field('p_UV_12')) # 'p_UV_13', 'p_UV_17' - other integration limits
p = interp1d(Red_z,rho_UV,kind='linear',bounds_error=None) # p is rho UV(z)

############################################################
# Integration limits and conversion z to time
#
ts = np.linspace(0.25,13.4,10000000)  # Produce time grid
zs= t_to_z(ts)
t0 = 0.0
###############################################################


#### START PLOT

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

#################################################################
# LAEs:
xi_ion = 10.**(25.6) #production_efficiency [Hz erg**(-1.)]
f_esc0 = 0.13

###### Compute how much UV luminosity density for LAEs ######
z_frac =[0,2,4,6,8,10,12,16]
fract_UV = [0.02,0.05,0.2,0.5,0.6,0.7,0.7,0.7]  # LAEs representing ~2-5% at z~0-2
leftovers = 1.0-array(fract_UV) # Track the remaining luminosity density
fract_for_fit = p(array(z_frac))*array(fract_UV)
leftover_lum = p(array(z_frac))*array(leftovers)

fract_UUV = interp1d(z_frac,fract_for_fit,kind='linear',bounds_error=None) # p is rho UV(z)
p = fract_UUV # Use these fractional rho_UVs as the luminosity densities

Q = odeint(dQ_dt,t0,ts)  # SOLVE NOW
Q_LAE=Q
Q[Q>1.0] = 1.0 # boundary condition up
Q[Q<0.0] = 0.0 # boundary condition down

LAEs = plot(zs,Q,'b-',linewidth=2,) # PLOT

#################################################################
# LBGs: non LA-emitting
xi_ion = 10.**(25.4) #production_efficiency [Hz erg**(-1.)]
f_esc0 = 0.03 #limits from stacking dominant non-LAEs at z~3

p =  interp1d(z_frac,leftover_lum,kind='linear',bounds_error=None) # p is rho UV(z)

Q = odeint(dQ_dt,t0,ts)  # SOLVE NOW
Q_LBG=Q
Q[Q>1.0] = 1.0 # boundary condition up
Q[Q<0.0] = 0.0 # boundary condition down

LBGS = plot(zs,Q,'r-',linewidth=3) # PLOT


#### TOTAL ###################
#
Q_TOT = Q_LAE+Q_LBG
Q_TOT[Q_TOT>1.0]=1.0
TOTAL = plot(zs,Q_TOT,'g--',linewidth=4,zorder=-1)


# LEGEND
I1=legend((TOTAL[0],LAEs[0],LBGS[0]),(r'Total LBGs',r'LAEs',r'Non-LAEs'), shadow = True, loc=3,numpoints=1)
ltext = gca().get_legend().get_texts()
setp(ltext[0], fontsize = 24, color='0.1')
setp(ltext[1], fontsize = 24, color='k')
#setp(ltext[2], fontsize = 15, color = '0.1')

frame=I1.get_frame()
frame.set_linewidth(0)
frame.set_visible(False)

# Define the axis
axis([0,12.1,-0.05,1.05])

# Label axis
xlabel(r'Redshift ($z$)', {'color':'k', 'fontsize': 23})
ylabel(r'Ionised Hydrogen Fraction (\%)', {'color':'k', 'fontsize': 23})

xticks((0,2,4,6,8,10,12),(0,2,4,6,8,10,12), color='k', size=22)
yticks((0.0,0.2,0.4,0.6,0.8,1.0),(0,20,40,60,80,100), color='k', size=22)

savefig('ionised_fraction.ps')
show()

