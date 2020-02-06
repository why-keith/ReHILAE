import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math

#TIME VALUES
startTime=0.01
finishTime=14
timeStep=0.01

# FUNCTION TO GENERATE Z AND T ARRAYS
def redshift(startTime, finishTime, timeStep): 
    
    #global t
    t = np.arange(startTime, finishTime, timeStep)
    Z = []
    
    for i in t:
        Z.append(mod.z(i))
    
    return Z,t

def oldRedshift(startTime, finishTime, timeStep):
    
    #global t
    t = np.arange(startTime, finishTime, timeStep)
    Z_old = []
    
    for i in t:
        Z_old.append(mod.z_old(i))
    
    return Z_old, t


#QUICK PLOT
    """
Z = redshift(startTime, finishTime, timeStep)
plotter.plot(t,Z,"A plot of redshift against time","z","t")
"""