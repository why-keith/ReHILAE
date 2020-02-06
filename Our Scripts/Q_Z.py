import numpy as np
from scipy.integrate import odeint
import MPhys_model
import matplotlib.pyplot as plt

ts = np.linspace(0.051,14,10000000) # time in Gyr
zs= ((((28./(ts))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift

def Q_Hii_dot(Q,t):
    z= ((((28./(t))-1.)**(1./2.)-1.))
    dQ_dt = MPhys_model.Q_Hii_dot(z,Q)
    return dQ_dt

Q = odeint(Q_Hii_dot, 0., ts)
Q[Q>1.0] = 1.0 # 100% HII
Q[Q<0.0] = 0.0 # 100% HI

plt.plot(zs,Q)
plt.xlabel('Redshift')
plt.ylabel('Fraction of Ionised Hydrogen')
plt.title('Graph to show the evolution of ionised hydrgoen against redshift')
plt.show()