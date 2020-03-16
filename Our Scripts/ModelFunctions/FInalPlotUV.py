import matplotlib.pyplot as plt
import numpy as np
import math
import random_array_generator as rag
import pandas as pd
import UV_Model as main
import time

start=time.time()

ts = np.linspace(0.051,14,100000) # time in Gyr
zs= ((((28./(ts))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift

C1 = [0.0013679194302549992]
C1_error = [0.0009034769902604823]

C2 = [-0.05592731916980156]
C2_error = [0.024847250273459104]

C3 = [0.29377928649363305]
C3_error = [0.21010382042922288]

C4 = [26.095271603086044]
C4_error = [0.5405447975405038]
print("Generating arrays...")
C1_P_UV = rag.random_Arrays(len(C1),C1,C1_error,C1_error)
C2_P_UV = rag.random_Arrays(len(C2),C2,C2_error,C2_error)
C3_P_UV = rag.random_Arrays(len(C3),C3,C3_error,C3_error)
C4_P_UV = rag.random_Arrays(len(C4),C4,C4_error,C4_error)

data = []
for i,j,k,l in zip(C1_P_UV, C2_P_UV, C3_P_UV, C4_P_UV):
    arguements = (i[0], j[0], k[0], l[0])
    data.append((main.main(ts,arguements)))

#print('simualtion finished')
selectedData=[]
for result in data:
    anonmalies = result[:13]
    if any([q[0]==1 for q in anonmalies]):
        continue
    else:
        selectedData.append(result)
#print('data filtered')
#print(len(data))
#print(len(selectedData))
median, median_lower_percentile, median_upper_percentile = rag.median_y_values(len(selectedData[0]),selectedData)
plt.figure()
for result in selectedData:
    plt.plot(zs,result)
plt.xlabel("Redshift (z)")
plt.ylabel("Fractions of Ionised Hydrogen")





plt.figure('Ionised_Hydrogen_UV')
plt.xlabel("Redshift (z)")
plt.ylabel("Fractions of Ionised Hydrogen")
plt.plot(zs,median, color = "black", label="LAE")

plt.fill_between(zs,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "steelblue", edgecolor = "black", linewidth = 1.2)
plt.fill_betweenx(median,6,10, color = "lightgrey", alpha = 0.3, edgecolor = "black", linewidth = 5)

plt.show()
print("Time elapsed = {}s".format(round(time.time()-start,2)))