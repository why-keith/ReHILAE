import numpy as np
from scipy.integrate import odeint
import LyCmodel as model
import matplotlib.pyplot as plt

def dQ_dt(Q,t):
    z= ((((29./(t))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift
    dQ_dt = model.Q_Hii_dot(z,Q)
    return dQ_dt

ts = np.linspace(0.051,14,10000000) # time in Gyr
zs= ((((29./(ts))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift

Q = odeint(dQ_dt, 0., ts) 
Q[Q>1.0] = 1.0 # 100% HII boundary condition up
Q[Q<0.0] = 0.0 # 100% HII boundary condition down

plt.figure('Ioinsed Hydrogen Fraction')
plt.plot(zs, Q)
plt.title('A Plot of the Fraction of Ionised Hydrogen Against Redshift')
plt.xlabel('Redshift ($z$)')
plt.ylabel('Fraction of Ionised Hydrogen ($Q_{II}$)')
plt.show()