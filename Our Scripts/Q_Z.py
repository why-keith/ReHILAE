import numpy as np
from scipy.integrate import odeint
import MPhys_model as mod
import matplotlib.pyplot as plt
from Redshift import redshift as z

#Z VALUES
startZ = mod.startZ
finishZ = mod.finishZ
zStep = mod.zStep

#GENERATE Z AND t ARRAYS
Z,t = z(startZ, finishZ, zStep)

#RETURNS dQ_dt 
def Q_Hii_dot(Q,t):
    dQ_dt = mod.Q_Hii_dot(mod.z(t),Q)
    return dQ_dt

#GENERATE Q ARRAY
Q = odeint(Q_Hii_dot, mod.Q_Hii_zero, t)
Q[Q>1.0] = 1.0 # 100% HII
Q[Q<0.0] = 0.0 # 100% HI

#PLOT Z AGAINST Q
plt.plot(Z,Q)
plt.xlabel('Redshift')
plt.ylabel('Fraction of Ionised Hydrogen')
plt.title('Graph to show the evolution of ionised hydrgoen against redshift')
plt.show()