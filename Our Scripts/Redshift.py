import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math

#TIME VALUES
startTime = mod.startTime
finishTime = mod.finishTime
timeStep = mod.timeStep

#GENERATES Z AND T ARRAYS
def redshift(startTime, finishTime, timeStep, alt = False): 
    
    t = np.arange(startTime, finishTime, timeStep)
    Z = []
    
    for i in t:
        Z.append(mod.z(i, alt))
    
    return Z,t
