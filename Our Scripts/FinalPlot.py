import matplotlib.pyplot as plt
import numpy as np
import math
import random_array_generator as rag
import pandas as pd
import Main_LAE
from Redshift import redshift

#TIME CONDITIONS
startT = Main_LAE.startT
finishT = Main_LAE.finishT
TStep = Main_LAE.TStep



"""
Power Law
"""


######################################################
#GENERATE P_L_Lya ITERATIONS
P1_1 = [1.13869]
P1_1_error = [0.18238]

P1_2 = [39.18853]
P1_2_error = [0.10799]

 
P1_1_P_L_Lya = rag.random_Arrays(len(P1_1),P1_1,P1_1_error,P1_1_error)
P1_2_P_L_Lya = rag.random_Arrays(len(P1_2),P1_2 ,P1_2_error,P1_2_error)
"""
P_L_Lya=[]
for i in range(len(A_P_L_Lya)):
    temp = []
    for j in range(len(z)):
        temp.append((A_P_L_Lya[i][0])*(z[j])*(z[j]) + (B_P_L_Lya[i][0])*z[j] + (C_P_L_Lya[i][0]))
    P_L_Lya.append(temp)
    
for i in P_L_Lya:
    plt.plot(z,i)
plt.show()
"""
####################################################### 
#GENERATE f_esc ITERATIONS

f1 = [0.00064]
f1_error = [0.00013]

f2 = [0.00941]
f2_error = [0.00364]
 
f1_f_esc = rag.random_Arrays(len(f1),f1,f1_error,f1_error)
f2_f_esc = rag.random_Arrays(len(f2),f2,f2_error,f2_error)
"""
f_esc=[]
for i in range(len(A_f_esc)):
    temp = []
    for j in range(len(z)):
        temp.append((A_f_esc[i][0])*(z[j]) + (B_f_esc[i][0]))
    f_esc.append(temp)
    
for i in f_esc:
    plt.plot(z,i)
plt.show()
"""
"""
#######################################################
data = []
z,t = redshift(startT, finishT, TStep)
for i,j,k,l in zip(P1_1_P_L_Lya,P1_2_P_L_Lya,f1_f_esc,f2_f_esc):
    arguements = (i[0],j[0],k[0],l[0])
    data.append((Main_LAE.main_PL(arguements))) # TODO write this return to file and then plot after loop


plt.figure('Ionised_Hydrogen_10')

#######################################################
#PLOTS MEDIAN OF Q_Hii_dot AND SHADES PERCENTILES
median, median_lower_percentile, median_upper_percentile = rag.median_y_values(len(data[0]),data)
plt.xlabel("Redshift (z)")
plt.ylabel("Fractions of Ionised Hydrogen")
plt.plot(z,median, color = "midnightblue", label="Median_PL")
plt.legend()

plt.fill_between(z,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "steelblue", edgecolor = "black", linewidth = 1.2)
plt.fill_betweenx(median,6,10, color = "skyblue", alpha = 0.3, edgecolor = "black", linewidth = 1.2)
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

for i in range(len(temp2)):
    plt.plot(z,data[temp2[i]])

plt.xlabel('Redshift (z)')
plt.ylabel('Fraction of Ionised Hydrogen')
"""
   


"""
Cubic
"""
    
######################################################
#GENERATE P_L_Lya ITERATIONS
P1_1 = [-0.04764]
P1_1_error = [0.00254]

P1_2 = [0.45297]
P1_2_error = [0.03671]

P1_3 = [38.97568]
P1_3_error = [0.08654]


 
P1_1_P_L_Lya = rag.random_Arrays(len(P1_1),P1_1,P1_1_error,P1_1_error)
P1_2_P_L_Lya = rag.random_Arrays(len(P1_2),P1_2 ,P1_2_error,P1_2_error)
P1_3_P_L_Lya = rag.random_Arrays(len(P1_3),P1_3,P1_3_error,P1_3_error)

"""
P_L_Lya=[]
for i in range(len(A_P_L_Lya)):
    temp = []
    for j in range(len(z)):
        temp.append((A_P_L_Lya[i][0])*(z[j])*(z[j]) + (B_P_L_Lya[i][0])*z[j] + (C_P_L_Lya[i][0]))
    P_L_Lya.append(temp)
    
for i in P_L_Lya:
    plt.plot(z,i)
plt.show()
"""
####################################################### 
#GENERATE f_esc ITERATIONS

f1 = [0.00064]
f1_error = [0.00013]

f2 = [0.00941]
f2_error = [0.00364]
 
f1_f_esc = rag.random_Arrays(len(f1),f1,f1_error,f1_error)
f2_f_esc = rag.random_Arrays(len(f2),f2,f2_error,f2_error)
"""
f_esc=[]
for i in range(len(A_f_esc)):
    temp = []
    for j in range(len(z)):
        temp.append((A_f_esc[i][0])*(z[j]) + (B_f_esc[i][0]))
    f_esc.append(temp)
    
for i in f_esc:
    plt.plot(z,i)
plt.show()
"""
#######################################################
data = []
z,t = redshift(startT, finishT, TStep)
for i,j,k,l,m in zip(P1_1_P_L_Lya,P1_2_P_L_Lya,P1_3_P_L_Lya,f1_f_esc,f2_f_esc):
    arguements = (i[0],j[0],k[0],l[0], m[0])
    data.append((Main_LAE.main_Q(arguements))) # TODO write this return to file and then plot after loop


plt.figure('Ionised_Hydrogen_10')

#######################################################
#PLOTS MEDIAN OF Q_Hii_dot AND SHADES PERCENTILES
median, median_lower_percentile, median_upper_percentile = rag.median_y_values(len(data[0]),data)
plt.xlabel("Redshift (z)")
plt.ylabel("Fractions of Ionised Hydrogen")
plt.plot(z,median, color = "black", label="Median_cubic")
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
plt.show()

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

for i in range(len(temp2)):
    plt.plot(z,data[temp2[i]])

plt.xlabel('Redshift (z)')
plt.ylabel('Fraction of Ionised Hydrogen')
plt.show()          
    
    
    
    