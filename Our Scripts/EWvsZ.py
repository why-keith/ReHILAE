"""
Generates a plot of the equivalent width against redshift
"""
import matplotlib as plt
import numpy as np
import MPhys_model as mod
import plotter
import math
from Redshift import redshift as z
import plot_saver as save
import random_array_generator as rag
import pandas as pd

if __name__=="__main__":
    path=None
else:
    path=save.folder

#TIME CONDITIONS
startT = mod.startT
finishT = mod.finishT
TStep = mod.TStep

#GENERATE Z AND t ARRAYS
z,t = z(startT, finishT, TStep)
z = np.asarray(z)

#GENERATES EW ARRAY
EW=[mod.EW(i) for i in t]

#PLOTS z AGAINST EW
plotter.plot(z,EW,"A plot of the equivalent width against redshift","Redshift","EW",path)

#DUMMY ERRORS
error_down_y = []
error_up_y = []

for i in range(len(z)):
        error_down_y.append(EW[i]*0.05)
        error_up_y.append(EW[i]*0.05)


#GENERATES AND SAVES RANDOM ARRAYS 
y = rag.random_Arrays(z,EW,error_down_y,error_up_y)
x = z
Data = [x,y]
df =  pd.DataFrame([Data], columns = ['Z', 'EW'])
df.to_pickle("EW_data.csv")
df = pd.read_pickle("EW_data.csv")
print(df)
