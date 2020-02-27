"""
Generates a plot of the evolution of ionised hydrgoen against redshift
"""
import numpy as np
from scipy.integrate import odeint
import AllModelsUV_Framework as mod
import matplotlib.pyplot as plt

#Z VALUES
ts = np.linspace(0.051,14,100000)
zs= ((((29./(ts))-1.)**(1./2.)-1.))

#RETURNS dQ_dt 
def dQ_dt_UV(Q,t):
    dQ_dt = mod.Q_Hii_dot_UV(mod.z(t),Q)
    return dQ_dt

def dQ_dt_LBG(Q,t):
    dQ_dt = mod.Q_Hii_dot_LBG(mod.z(t),Q)
    return dQ_dt

def dQ_dt_LAE(Q,t):
    dQ_dt = mod.Q_Hii_dot_Lya(mod.z(t),Q)
    return dQ_dt

#GENERATE Q ARRAY
Q_UV = odeint(dQ_dt_UV, mod.Q_Hii_zero, ts)
Q_UV[Q_UV>1.0] = 1.0 # 100% HII
Q_UV[Q_UV<0.0] = 0.0 # 100% HI

Q_LBG = odeint(dQ_dt_LBG, mod.Q_Hii_zero, ts)
Q_LBG[Q_LBG>1.0] = 1.0 # 100% HII
Q_LBG[Q_LBG<0.0] = 0.0 # 100% HI

Q_Lyc = odeint(dQ_dt_LAE, mod.Q_Hii_zero, ts)
Q_Lyc[Q_Lyc>1.0] = 1.0 # 100% HII
Q_Lyc[Q_Lyc<0.0] = 0.0 # 100% HI

#PLOT Z AGAINST Q
plt.figure('UV Framework')
plt.plot(zs,Q_UV, color='green', label='Robertson Model')
plt.plot(zs,Q_LBG, '-', color='red', label='LBGs')
plt.plot(zs,Q_Lyc, '-', color='blue', label='LAEs')
plt.xlabel('Redshift ($z$)')
plt.ylabel('Fraction of Ioinsed Hydrogen')
plt.title('Modelling the Epoch of Reionisation')
plt.legend()
plt.show()
