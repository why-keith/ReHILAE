"""
Produces arrays of redshift and time based on a start 
and finish redshift and an incremental step size
"""
import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math
import random_array_generator as rag
import pandas as pd


#Z VALUES
#startZ = mod.startZ
#finishZ = mod.finishZ
#zStep = mod.zStep

#TIME CONDITIONS
startT = mod.startT
finishT = mod.finishT
TStep = mod.TStep

#GENERATES Z, Z_old, AND T ARRAYS
def redshift_old(startZ, finishZ, zStep, alt = False): 
    
    Z = np.arange(startZ, finishZ, zStep)
    T = []
    
    for i in Z:
        T.append(mod.t(i, alt))
    
    return Z,T

def redshift(startT, finishT, TStep, alt = False):
    
    T = np.arange(startT, finishT, TStep)
    Z=mod.z(T, alt)
    return Z,T