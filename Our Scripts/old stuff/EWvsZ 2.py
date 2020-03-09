"""
Generates a plot of the equivalent width against redshift
"""
import matplotlib.pyplot as plt
import pylab
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
EW=[mod.EW(i) for i in z]


x = pylab.array([2.5,2.8,2.9,3.1, 3.3, 3.7, 4.1, 4.5, 4.8, 5.0, 5.3])
y = pylab.array([117.134349876, 122.113769531, 172.679885864, 143.5936203, 149.371124268, 96.3648490905, 189.002449036, 181.127731323, 95.0448436864, 125.935153962, 143.343048096])

#PLOTS z AGAINST EW
plt.scatter(x, y)
plt.plot(z, EW )
plt.title("A plot of the equivalent width against redshift")
plt.xlabel("Redshift")
plt.ylabel("EW [Å]")
plt.show()


#plotter.plot(z,EW,"A plot of the equivalent width against redshift","Redshift","EW [Å]",path) 

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