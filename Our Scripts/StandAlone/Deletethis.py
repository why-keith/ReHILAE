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

#
### Main properties/constants
#
xi_ion = 10.**(25.6) #production_efficiency [Hz erg**(-1.)]
f_esc0 = 0.2
alpha = 1.17 # this is for f_esc dependence on redshift
x1 = array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4,17])
y1 = array([10.**(26.52), 10.**(26.30), 10.**(26.10), 10.**(25.98), 10.**(25.67), 10.**(25.32),10.**(23.12) ])
p1 = polyfit(x1,y1,3.)
p = poly1d(p1)


Omega_b = 0.04 #baryon_density
h = 0.7
X_p = 0.75 #hydrogen_mass_fraction
T = 20000. #IGM_temperature [K]
C = 5  #clumping_factor
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



### DEFINE SIZE

#figsize=(12,12)
#f, (ax1) = plt.subplots(1, sharey=True, figsize=(7.3,6.6))


########################################
# Plot end of re-ionisation

za = [6.0,6.0000001]
zp = [-0.5,5.0]

Reionisation = plot(za,zp,'k--')

## PLOT REIONISATION EPOCH
#
za = [8.8,8.80001]
zp = [-0.5,5.0]

Reionisation = plot(za,zp,linestyle='-',color ='0.40', linewidth=8, zorder=-8)

#Reionisation = plot(array(za)-1.1,zp,linestyle='-.',color ='k', linewidth=4, zorder=-8) lower limit inst reion
#Reionisation = plot(array(za)+1.2,zp,linestyle='-.',color ='k', linewidth=4, zorder=-8) upper limit inst reion

Reionisation = plot(array(za)+0.05,zp,linestyle='-',color ='0.90', linewidth=100, zorder=-10)

Reionisation = plot(array(za)+0.05,zp,linestyle='-',color ='0.97', linewidth=200, zorder=-20)

#######################################
# Integration limits and conversion z to time
#
ts = np.linspace(0.051,14,10000000)
zs= ((((28./(ts))-1.)**(1./2.)-1.))
t0 = 0.0
#

#################################################################
# FOR LBGs:
xi_ion = 10.**(25.4) #production_efficiency [Hz erg**(-1.)]
f_esc0 = 0.07
x1 = array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4,17])
y1 = array([10.**(26.52), 10.**(26.30), 10.**(26.10), 10.**(25.98), 10.**(25.67), 10.**(25.32),10.**(23.12) ])
p1 = polyfit(x1,y1,3.)
p = poly1d(p1)

Q = odeint(dQ_dt,t0,ts)  # SOLVE NOW
Q[Q>1.0] = 1.0
Q[Q<0.0] =0.0
LBGS = plot(zs,Q,'r-',linewidth=4,zorder = -2) # PLOT

Q_LBG = Q

###############################################################
# FOR LAEs:
xi_ion = 10.**(25.65) #production_efficiency [Hz erg**(-1.)]
f_esc0 = 0.20
x1 = array([3.8, 4.9, 5.9, 6.8, 7.9, 10.4,17])
y1 = array([10.**(26.52), 10.**(26.30), 10.**(26.10), 10.**(25.98), 10.**(25.67), 10.**(25.32),10.**(23.12) ])

factors = array([0.40, 0.30, 0.5, 0.8, 1.0, 1.0,1.0 ]) # Fraction of UV luminosity density recovered by Lyman-alpha emitters
y1 = y1 * factors

p1 = polyfit(x1,y1,3.)
p = poly1d(p1)

Q = odeint(dQ_dt,t0,ts)
Q[Q>1.0] = 1.0
Q[Q<0.0] =0.0
Q[zs<3] = 1.0

Q_LAE = Q
LAEs = plot(zs,Q,'b--',linewidth=6,zorder=1)

#### TOTAL
#
#
Q_TOT = Q_LAE+Q_LBG
Q_TOT[Q_TOT>1.0]=1.0
TOTAL = plot(zs,Q_TOT,'g-',linewidth=3,zorder=-1)



# LEGEND
I1=legend((TOTAL[0],LBGS[0],LAEs[0]),(r'Total',r'LBGs',r'LAEs'), shadow = True, loc=6,numpoints=1)
ltext = gca().get_legend().get_texts()
setp(ltext[0], fontsize = 24, color='0.1')
setp(ltext[1], fontsize = 24, color='k')
#setp(ltext[2], fontsize = 15, color = '0.1')

frame=I1.get_frame()
frame.set_linewidth(0)
frame.set_visible(False)


# Define the axis
axis([0,11.1,-0.05,1.05])

# Label axis
xlabel(r'Redshift ($z$)', {'color':'k', 'fontsize': 23})
ylabel(r'Ionised Hydrogen Fraction (\%)', {'color':'k', 'fontsize': 23})

xticks((0,2,4,6,8,10),(0,2,4,6,8,10), color='k', size=22)
yticks((0.0,0.2,0.4,0.6,0.8,1.0),(0,20,40,60,80,100), color='k', size=22)

show()