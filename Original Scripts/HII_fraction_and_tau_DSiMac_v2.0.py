##############################################
#
# Solving Re-ionisation
# By DS 2017, updated by DS 2018
#
##############################################
from pylab import *
from scipy.integrate import odeint
import matplotlib
#matplotlib.rc('text', usetex=True)   # Commented out for iMAC
from math import *
from scipy.interpolate import interp1d
import pyfits as py
from scipy.integrate import quadpack
from scipy.integrate import quadrature
from scipy import special
import os

# Quick function to open fits catalogue and deliver it back
def fits_open(file_in):
    Open_C=pyfits.open(file_in)
    Fits_catalogue=Open_C[1].data
    return Fits_catalogue


#### LOAD SC4K RESULTS
def Load_SC4K(sample='SSC4K'):

    CAT='../Santos2017_SC4K/Analysis_Global_Lya_LF/MCMC_Results_11122017_vfinal2.fits'
    
    # Open the catalogue
    data=fits_open(CAT)
    # FIELDS with relevant information
    SUB_SET = data.field('Sub_set_condition') # Name of sub-set
    Alpha = data.field('Alpha') # Name of sub-set
    Upper_lum_fit = data.field('Fit_upp_luminosity')
    
    #LOAD ALL info
    L_dens, Ldens_up, Ldens_down, SFRD_, SFRD_up, SFRD_down, Ldens_PL, Ldens_PL_up, Ldens_PL_down = data.field('Ldens'), data.field('Ldens_error_up'), data.field('Ldens_error_down'), data.field('Ldens')*Dens_to_SFRD, data.field('Ldens_error_up')*Dens_to_SFRD, data.field('Ldens_error_down')*Dens_to_SFRD, data.field('Ldens_PL2'), data.field('Ldens_error_up_PL2'), data.field('Ldens_error_down_PL2')
    # SETUP ERRORS
    ####
    Ldens_up_log   =  np.log10(L_dens+Ldens_up)-np.log10(L_dens)
    Ldens_down_log = np.log10(L_dens)-np.log10(L_dens-Ldens_down)
    #
    Ldens_PL_up_log =  np.log10(Ldens_PL+Ldens_PL_up)-np.log10(Ldens_PL)
    Ldens_PL_down_log =  np.log10(Ldens_PL)-np.log10(Ldens_PL-Ldens_PL_down)
    
    # GET REDSHIFTS
    Redshift_c, Delta_redshift = (data.field('Redshift_min')+data.field('Redshift_max'))/2.,(data.field('Redshift_max')-data.field('Redshift_min'))/2.
    # Copy to each filter
    MASK_each_filter = Alpha>100000.0
    MASK_each_fil_DRAKE= Alpha>100000.0
    # Now get the conditions properly:
    FILTERS=['NB392','IA427','IA464','IA484','IA505','IA527','IA574','IA624','IA679','IA709','IA738','IA767','IA827']
    INTEGRAL_con = [43.3,43.3,43.3,43.3,43.3,43.3,44.5,44.5,44.5,44.5,44.5,44.5,44.5]
    #
    for i in range(len(FILTERS)):
        MASK_each_filter = MASK_each_filter+ ((SUB_SET == '%s'%FILTERS[i])*ALPHA_fixed*(Upper_lum_fit==INTEGRAL_con[i]))
        MASK_each_fil_DRAKE = MASK_each_fil_DRAKE+ ((SUB_SET == '%s_with_All_lit'%FILTERS[i])*(Upper_lum_fit==INTEGRAL_con[i]))
    
    # NOW for redshift bins
    FILTERS = ['z22_ALL','z25_ALL', 'z31_ALL', 'z39_ALL', 'z47_ALL', 'z54_ALL']
    INTEGRAL_con = [43.3,43.3,43.3,44.5,44.5,44.5]
    MASK_stacks = Alpha>100000.0
    MASK_stacks_DRAKE= Alpha>100000.0
    #
    for i in range(len(FILTERS)):
        MASK_stacks = MASK_stacks+ ((SUB_SET == '%s'%FILTERS[i])*ALPHA_fixed*(Upper_lum_fit==INTEGRAL_con[i]))
        MASK_stacks_DRAKE = MASK_stacks_DRAKE+ ((SUB_SET == '%s_with_All_lit'%FILTERS[i])*(Upper_lum_fit==INTEGRAL_con[i]))
        
    # Now for full stack
    # NOW for redshift bins
    FILTERS = ['FULL_SC4K','FULL_SC4K']
    INTEGRAL_con = [43.3,43.5]
    FULL_SC4K = Alpha>100000.0
    FULL_SSC4K = Alpha>100000.0
    #
    for i in range(len(FILTERS)):
        FULL_SC4K = FULL_SC4K+ ((SUB_SET == '%s'%FILTERS[i])*ALPHA_fixed*(Upper_lum_fit==INTEGRAL_con[i]))
        FULL_SSC4K = FULL_SSC4K+ ((SUB_SET == '%s_with_All_lit'%FILTERS[i])*(Upper_lum_fit==INTEGRAL_con[i]))
    
    
    if sample=='SSC4K':
        return Redshift_c[MASK_each_fil_DRAKE],np.log10(L_dens[MASK_each_fil_DRAKE]), Ldens_down_log[MASK_each_fil_DRAKE],Ldens_up_log[MASK_each_fil_DRAKE]
    


############################################
# Return condicence levels
def plot_confidence_Lya(x_in,ysample,colour='g'):
    lower = percentile(ysample,16,axis=0)  # Equivalent to 1 sigma
    upper = percentile(ysample,84,axis=0) # Equivalent to 1 sigma
    lower2 = percentile(ysample,100-97.72,axis=0)  #, Equivalent to 2 sigma
    upper2 = percentile(ysample,97.72,axis=0) # Equivalent to 2 sigma
    lower3 = percentile(ysample,0.15,axis=0)  #, Equivalent to 3 sigma
    upper3 = percentile(ysample,99.85,axis=0) # Equivalent to 3 sigma
    return x_in,lower,upper, lower2, upper2, lower3, upper3

####################################
# Gaussian Function
def gaussian(x,x0,sigma):
    return np.exp(-(x-x0)**2/(2*sigma**2))

####################################
# Double gaussian, asymmetric
# phi_dummy generates the space of values to pick from
# phi_probs generates the PDF for a 1sigma error around each phi
def double_normal(phi,err_down,err_up,size):
   phi_dummy = np.linspace(phi-5.*err_down,phi+5.*err_up,100000) #generate list from -5 to +5 error
   phi_probs = np.append(gaussian(phi_dummy[phi_dummy<phi],phi,err_down),gaussian(phi_dummy[phi_dummy>=phi],phi,err_up))
   return np.random.choice(phi_dummy,size,p=phi_probs/np.sum(phi_probs))
   
##################################
def dQ_dt (Q,t): #HII_fraction

    z = t_to_z(t)
    trec = 1./(C*alpha_B*n_e*(1.+z)**(3.)) # s 
    nion = f_esc0*xi_ion*p(z)/ (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3
    return ((nion)/n_H-Q/trec)*3.1536e+16 # Conversion to time in seconds
    

##################################
def dQ_dt_Lya (Q,t): #HII_fraction

    z = t_to_z(t)
    trec = 1./(C*alpha_B*n_e*(1.+z)**(3.)) # s 
    cHa_K98 = 1.36*pow(10.,-12.)
    
    #Q_ionLya =  LLya/(cHa_K98*(1-f_esc0)(0.042*EW0)) # Equation 5 Sobral & Matthee 2018
    nion = f_esc0 * (pLya(z)/(cHa_K98*(1.-f_esc0)*(0.042*EW0(z))))/ (2.938e+73) # From Sobral & Matthee 2018
    #nion = f_esc0*xi_ion*p(z)/ (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3
    return ((nion)/n_H-Q/trec)*3.1536e+16 # Conversion to time in seconds


##################################
def dQ_dt_direct (Q,t): #HII_fraction

    z = t_to_z(t)
    trec = 1./(C*alpha_B*n_e*(1.+z)**(3.)) # s 
    cHa_K98 = 1.36*pow(10.,-12.)
    
    nion = pLyC(z)/ (2.938e+73) # (2.938e+73) converts from Mpc^-3 to cm^-3  - full units s^-1 Mpc^-3
    return ((nion)/n_H-Q/trec)*3.1536e+16 # Conversion to time in seconds
    
    
##################################
def tau (z):

    return sigma_T*n_e*(1.+z)*(1.+z)*Q(z)*c/(H)/z   # Units:
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
c = 299792.458  # Speed of light km/s

#### For getting Tau (z)
sigma_T = 6.65*pow(10,-25) # Thomson scattering cross section   cm^2  - from Mayer et al. 
n_e = (1.+Y_p/(4.*X_p))*n_H   # comoving number density of electrons (assuming singly ionized He)
H = h*100. *(3.24077929*pow(10,-25.)) # Hubble parameter == hubble constant conversion per cm from per Mpc

#######################################################################################

##################################################################
# GET redshift - time conversion
######
# F737 Cosmology
F737 = 'Redshift_time_F737_cosmology.fits'
# PLanck cosmology
FPlanck = 'Redshift_time_FPlanck_cosmology.fits'

# Which Cosmology to use:
USING_NOW = FPlanck

#LOAD Redshift-time dependence
A=py.open(USING_NOW) 
NOW=A[1].data
#LOAD GENERAL PROPERTIES - what will be plotted
Redshift = NOW.field('Redshift') 
Time =  NOW.field('Time_Gyr') 
# time to redshift relation
t_to_z = interp1d(Time,Redshift,kind='linear',bounds_error=None)

############################################################
# Integration limits and conversion z to time
#
ts = np.linspace(0.25,13.4,50000)  # Produce time grid  - 50000 steps seems to be a good minimum for convergence
zs= t_to_z(ts)
t0 = 0.0
###############################################################



#################################

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
#Reionisation = plot(array(za)+0.05,zp,linestyle='-',color ='0.90', linewidth=100, zorder=-10)
#Reionisation = plot(array(za)+0.05,zp,linestyle='-',color ='0.97', linewidth=200, zorder=-20)

Reionisation = ax1.plot(za,zp,linestyle='-',color ='0.60', linewidth=8, zorder=-18)
#Reionisation = ax1.plot(array(za)+0.05,zp,linestyle='-',color ='0.94', linewidth=100, zorder=-10)
#Reionisation = ax1.plot(array(za)+0.05,zp,linestyle='-',color ='0.99', linewidth=200, zorder=-20)

############################



########################################
# UV luminosity densities
B=py.open('RHO_UV.fits') 
NOWB=B[1].data
Red_z = NOWB.field('Redshift') 
rho_UV =  10.**(NOWB.field('p_UV_12')) # 'p_UV_12  p_UV_13', 'p_UV_17' - other integration limits
rho_UV_log = (NOWB.field('p_UV_12')) # 'p_UV_12  p_UV_13', 'p_UV_17' - other integration limits
p = interp1d(Red_z,rho_UV,kind='linear',bounds_error=None) # p is rho UV(z)

# LOAD Lyman-alpha luminosity density, from S-SC4K
z_Lya, logLumLya, err_downLya, err_upLya=Load_SC4K(sample='SSC4K')
# Convert to list so can be appended
z_Lya, logLumLya, err_downLya, err_upLya= list(z_Lya), list(logLumLya), list(err_downLya), list(err_upLya)

##############################################################
# Add z=6.6, z=7.0, z=7.3 from Matthee+, Zheng+, Konno+
z_others = [6.6,7.0,7.3]
pLya_others = [39.4,38.92,38.55]
err_others = [0.3,0.4,0.3]
# Best linear fit
Bfit = [-1.21351351, 47.41081081]
# Best linear fit is given by: A=-1.21351351,  B=47.41081081   in log10
# Add higher redshift data - underestimation likely due to reionisaiton effects
for i in range(z_others):
    z_Lya.append(z_others[i])
    logLumLya.append(pLya_others[i])
    err_downLya.append(err_others[i])
    err_upLya.append(err_others[i])
# Add higher redshift extrapolated data
for i in [7.5,8.0,8.5,9.0,10.0,10.5,11.0,12.,13.0,14.,15.0]:
    z_Lya.append(i)
    logLumLya.append(i*Bfit[0]+Bfit[1])
    err_downLya.append(0.4)
    err_upLya.append(0.6)    
# Convert back to arrays
z_Lya, logLumLya, err_downLya, err_upLya= array(z_Lya), array(logLumLya), array(err_downLya), array(err_upLya)
# Fit higher redshift to extrapolate:


Q_storing, Tau_storing = [], []

iterations = 100

xi_ions = [25.2,0.2,0.6] # central value, 1 sigma down, 1 sigma up
xi_ion_lya = [25.5,0.2,0.4]
f_esc0_lya = [0.15,0.05,0.05]

f_esc0s = [0.15,0.075,0.05]
pUV_uncert = [0.1,0.05]
##### TRY MCMC quick

for i in range(iterations):#### Draw from normal distribution
            
            rho_lya_new = logLumLya-logLumLya # create new array with zeros
            for j in range(len(logLumLya)): # Go through each one
                rho_lya_new[j]=pow(10.,double_normal(logLumLya[j],err_downLya[j],err_upLya[j],1)[0])
            
            pLya(z) = interp1d(z_Lya,rho_lya_new,kind='linear',bounds_error=None) 
            # for now assume EW fixed
            EW0(z) = interp1d(z_Lya,rho_lya_new-rho_lya_new+80.0,kind='linear',bounds_error=None)
                
            rho_UV_new = rho_UV_log-rho_UV_log # create array with zeros
            for j in range(len(rho_UV_log)): # Go through each one
                rho_UV_new[j]=pow(10.,double_normal(rho_UV_log[j],pUV_uncert[0],pUV_uncert[1],1)[0])
            p = interp1d(Red_z,rho_UV_new,kind='linear',bounds_error=None)   
                
            # Pick an x_ion - in log10
            xi_ion=pow(10.,double_normal(xi_ions[0],xi_ions[1],xi_ions[2],1)[0]) ###ASYMETRIC
            # Now pick and f_esc LyC
            f_esc0 = double_normal(f_esc0s[0],f_esc0s[1],f_esc0s[2],1)[0]
            
            # Solve equation
            Q = odeint(dQ_dt_Lya,t0,ts)  # SOLVE NOW
            
            # SET boundary conditions
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
            Q_storing.append(Q(zs))
            
            # NOW TAU:
            
            z_plot, taus_plot = compute_tau_z(Q)
            
            Tau_storing.append(taus_plot)
            
            #plot(z_not,Q(z_not),'g-',linewidth=0.4,zorder=-1)

# FIND AND PLOT COUNTOURS FOR Fraction of Ionised Hydrogen
xi,lower,upper, lower2, upper2, lower3, upper3 = plot_confidence_Lya(zs,Q_storing,colour='#ff8080')
SHADE = fill_between(xi,lower,upper,lw=1,color='#0066ff',alpha=0.4,zorder = 90)
#SHADE = fill_between(xi,lower2,upper2,lw=1,color='#e6e6ff',alpha=0.8,zorder = 80)
#SHADE = fill_between(xi,lower3,upper3,lw=1,color='#cce6ff',alpha=0.4,zorder = -100)

# FIND AND PLOT CONTOUNRS FOR Tau - CMB optical depth to compare with Planck
xi,lower,upper, lower2, upper2, lower3, upper3 = plot_confidence_Lya(z_plot,Tau_storing,colour='#ff8080')
SHADE = ax1.fill_between(xi,lower,upper,lw=1,color='#0066ff',alpha=0.4,zorder = 90)




## TEST FOR LAEs 



# Solve equation
Q = odeint(dQ_dt_Lya,t0,ts)  # SOLVE NOW

# SET boundary conditions
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
            
            
               

########################################
# UV luminosity densities
B=py.open('RHO_UV.fits') 
NOWB=B[1].data
Red_z = NOWB.field('Redshift') 
rho_UV =  10.**(NOWB.field('p_UV_13')) # 'p_UV_12  p_UV_13', 'p_UV_17' - other integration limits
rho_UV_log = (NOWB.field('p_UV_13')) # 'p_UV_12  p_UV_13', 'p_UV_17' - other integration limits
p = interp1d(Red_z,rho_UV,kind='linear',bounds_error=None) # p is rho UV(z)



            
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


##################################################################################
### NOW TAU PLOT

# SOLVE TAU:
z_plot, taus_plot = compute_tau_z(Q)









TOTAL = plot(z_not,Q(z_not),'g-',linewidth=3,zorder=-1)
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

savefig('Q_HII_and_Tau_evolution.pdf')
show()
