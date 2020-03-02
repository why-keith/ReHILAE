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


######################################################
#GENERATE P_L_Lya ITERATIONS
a = [-0.05]
a_error = [0.01]

b = [0.44]
b_error = [0.09]

c = [38.99]
c_error = [0.15]
 
A_P_L_Lya = rag.random_Arrays(len(a),a,a_error,a_error)
B_P_L_Lya = rag.random_Arrays(len(b),b,b_error,b_error)
C_P_L_Lya = rag.random_Arrays(len(c),c,c_error,c_error)
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

a = [0.00064]
a_error = [0.00013]

b = [0.00941]
b_error = [0.00364]
 
A_f_esc = rag.random_Arrays(len(a),a,a_error,a_error)
B_f_esc = rag.random_Arrays(len(b),b,b_error,b_error)
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
for i,j,k,l,m in zip(A_P_L_Lya,B_P_L_Lya,C_P_L_Lya,A_f_esc,B_f_esc):
    arguements = (i[0],j[0],k[0],l[0],m[0])
    data.append((Main_LAE.main(arguements))) # TODO write this return to file and then plot after loop

#######################################################
#PLOTS MEDIAN OF Q_Hii_dot AND SHADES PERCENTILES
median, median_lower_percentile, median_upper_percentile = rag.median_y_values(len(data[0]),data)
plt.xlabel("Redshift (z)")
plt.ylabel("Fractions of Ionised Hydrogen")
plt.plot(z,median, color = "midnightblue", label="Median")
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
        
    
    
    
    
    
    