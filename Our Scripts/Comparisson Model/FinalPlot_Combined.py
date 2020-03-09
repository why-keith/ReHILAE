import matplotlib.pyplot as plt
import numpy as np
import math
import random_array_generator as rag
import pandas as pd
import Main_LAE_Powerlaw as main
import Main_UV

startT = main.startT
finishT = main.finishT
TStep = main.TStep

P1_1 = [1.13869]
P1_1_error = [0.18238]

P1_2 = [39.18853]
P1_2_error = [0.10709]

P2_1 = [-5.288]
P2_1_error = [0.83404]

P2_2 = [44.5388067]
P2_2_error = [0.757]


P1_1_P_L_Lya = rag.random_Arrays(len(P1_1),P1_1,P1_1_error,P1_1_error)
P1_2_P_L_Lya = rag.random_Arrays(len(P1_2),P1_2 ,P1_2_error,P1_2_error)
P2_1_P_L_Lya = rag.random_Arrays(len(P2_1),P2_1,P2_1_error,P2_1_error)
P2_2_P_L_Lya = rag.random_Arrays(len(P2_2),P2_2 ,P2_2_error,P2_2_error)

f1 = [0.00064]
f1_error = [0.00013]

f2 = [0.00941]
f2_error = [0.00364]
 
f1_f_esc = rag.random_Arrays(len(f1),f1,f1_error,f1_error)
f2_f_esc = rag.random_Arrays(len(f2),f2,f2_error,f2_error)

data = []
z,t = main.redshift(startT, finishT, TStep)
for i,j, k, l, m, n in zip(P1_1_P_L_Lya,P1_2_P_L_Lya, P2_1_P_L_Lya, P2_2_P_L_Lya, f1_f_esc, f2_f_esc):
    arguements = (i[0],j[0], k[0], l[0], m[0], n[0])
    data.append((main.main_PL_2(arguements))) # TODO write this return to file and then plot after loop


plt.figure('Ionised_Hydrogen_UV  10000')
plt.title("Fraction of ionised H with respect to redshift")

#######################################################
#PLOTS MEDIAN OF Q_Hii_dot AND SHADES PERCENTILES
median, median_lower_percentile, median_upper_percentile = rag.median_y_values(len(data[0]),data)
plt.xlabel("Redshift (z)")
plt.ylabel("Fractions of Ionised Hydrogen")
plt.plot(z,median, color = "black", label="LAE")


plt.fill_between(z,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "lightgrey", edgecolor = "black", linewidth = 1.2)
plt.fill_betweenx(median,6,10, color = "lightgrey", alpha = 0.3, edgecolor = "black", linewidth = 1.2)
#######################################################
#PLOTS ALL OF THE ITERATIONS
"""
plt.figure('Ionised Hydrogen  10000')
plt.title('Individual fits for each iteration' )
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
plt.title("selected fits for reionisation between z = 6 ~ 10")
for i in range(len(temp2)):
    plt.plot(z,data[temp2[i]])

plt.xlabel('Redshift (z)')
plt.ylabel('Fraction of Ionised Hydrogen')
"""
    
startT = Main_UV.startT
finishT = Main_UV.finishT
TStep = Main_UV.TStep

P1_1 = [-5.288]
P1_1_error = [0.83404]

P1_2 = [30.39894]
P1_2_error = [0.757]

P1_1_P_L_Lya = rag.random_Arrays(len(P1_1),P1_1,P1_1_error,P1_1_error)
P1_2_P_L_Lya = rag.random_Arrays(len(P1_2),P1_2 ,P1_2_error,P1_2_error)

data = []
z,t = Main_UV.redshift(startT, finishT, TStep)
for i,j in zip(P1_1_P_L_Lya,P1_2_P_L_Lya):
    arguements = (i[0],j[0])
    data.append((Main_UV.main_UV(arguements))) # TODO write this return to file and then plot after loop


plt.figure('Ionised_Hydrogen_UV  10000')
plt.title("Fraction of ionised H with respect to redshift")

#######################################################
#PLOTS MEDIAN OF Q_Hii_dot AND SHADES PERCENTILES
median, median_lower_percentile, median_upper_percentile = rag.median_y_values(len(data[0]),data)
plt.xlabel("Redshift (z)")
plt.ylabel("Fractions of Ionised Hydrogen")
plt.plot(z,median, color = "steelblue", label = "UV")
plt.legend()

plt.fill_between(z,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "lightblue", edgecolor = "black", linewidth = 1.2)
plt.fill_betweenx(median,6,10, color = "grey", alpha = 0.3, edgecolor = "black", linewidth = 1.2)
"""
#######################################################
#PLOTS ALL OF THE ITERATIONS
plt.figure('Ionised Hydrogen  10000')
plt.title('Individual fits for each iteration' )
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
plt.title("selected fits for reionisation between z = 6 ~ 10")
for i in range(len(temp2)):
    plt.plot(z,data[temp2[i]])

plt.xlabel('Redshift (z)')
plt.ylabel('Fraction of Ionised Hydrogen')
"""
plt.legend()
plt.show()          
    
    