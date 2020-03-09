import matplotlib.pyplot as plt
import numpy as np
import math
import random_array_generator as rag
import pandas as pd
import Main_LAE_Q as Main_LAE

#TIME CONDITIONS
startT = Main_LAE.startT
finishT = Main_LAE.finishT
TStep = Main_LAE.TStep

P1_1 = [-0.04764]
P1_1_error = [0.0128]

P1_2 = [0.45297]
P1_2_error = [0.0944]

P1_3 = [38.97568]
P1_3_error = [0.15874]
 
P1_1_P_L_Lya = rag.random_Arrays(len(P1_1),P1_1,P1_1_error,P1_1_error)
P1_2_P_L_Lya = rag.random_Arrays(len(P1_2),P1_2 ,P1_2_error,P1_2_error)
P1_3_P_L_Lya = rag.random_Arrays(len(P1_3),P1_3,P1_3_error,P1_3_error)

f1 = [0.00064]
f1_error = [0.00013]

f2 = [0.00941]
f2_error = [0.00364]
 
f1_f_esc = rag.random_Arrays(len(f1),f1,f1_error,f1_error)
f2_f_esc = rag.random_Arrays(len(f2),f2,f2_error,f2_error)

#######################################################
data = []
z,t = Main_LAE.redshift(startT, finishT, TStep)
for i,j,k,l, m in zip(P1_1_P_L_Lya,P1_2_P_L_Lya,P1_3_P_L_Lya, f1_f_esc,f2_f_esc):
    arguements = (i[0],j[0],k[0],l[0], m[0])
    data.append((Main_LAE.main_Q(arguements))) # TODO write this return to file and then plot after loop


plt.figure('Ionised_Hydrogen_LyC_Quad')

#######################################################
#PLOTS MEDIAN OF Q_Hii_dot AND SHADES PERCENTILES
median, median_lower_percentile, median_upper_percentile = rag.median_y_values(len(data[0]),data)
plt.xlabel("Redshift (z)")
plt.ylabel("Fractions of Ionised Hydrogen")
plt.plot(z,median, color = "black", label="Median_Quad")
plt.legend()

plt.fill_between(z,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "grey", edgecolor = "black", linewidth = 1.2)
plt.fill_betweenx(median,6,10, color = "lightgrey", alpha = 0.3, edgecolor = "black", linewidth = 1.2)
#######################################################
#PLOTS ALL OF THE ITERATIONS
plt.figure('Ionised Hydrogen')
plt.xlabel('Redshift (z)')
plt.ylabel('Fraction of Ionised Hydrogen')

for i in data:
    plt.plot(z,i)


######################################################
#FINDS ITERATIONS THAT ARE ALLOWED USING LAEs
temp1 = []
for i in range(len(data)):
    for j in range(len(data[i])):
        #print(z[j])
        if z[j] >=6 and z[j] <=10 and data[i][j-1] < 1 and data[i][j] ==1:
            temp1.append(i)
            break
  
temp2 = []
for i in range(len(temp1)):
    for j in range(len(data[i])):
        if z[j] >=6 and z[j] <=10 and data[temp1[i]][j] <= 0.1:
            temp2.append(temp1[i])
            break

print(len(temp1), len(temp2))

print("List of allowed iterations using LAEs only: " + str(temp2))
plt.figure('selected  10000')
for i in range(len(temp2)):
    plt.plot(z,data[temp2[i]])

plt.xlabel('Redshift (z)')
plt.ylabel('Fraction of Ionised Hydrogen')
plt.show()          
    
    
    
    