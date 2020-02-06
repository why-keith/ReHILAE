# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 11:05:40 2020

@author: philli20
"""

import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math
#import Redshift as z
from Redshift import redshift as z


#TIME VALUES
startTime=1
finishTime=6
timeStep = 0.01


#GENERATE Z ARRAY
#Z = z.redshift(startTime, finishTime, timeStep)
Z,t = z(startTime, finishTime, timeStep)

#GENERATE Z_alt ARRAY
#Z = z.redshift(startTime, finishTime, timeStep)
Z_alt,t = z(startTime, finishTime, timeStep, True)


#QUICK PLOT
plotter.plot(Z,Z_alt,"A plot of old redshift against new redshift","Redshift","Old Redshift")