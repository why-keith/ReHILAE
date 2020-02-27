import matplotlib.pyplot as plt
import numpy as np
import math
import random_array_generator as rag
import pandas as pd
import Main_LAE as main

#TIME CONDITIONS
startT = main.startT
finishT = main.finishT
TStep = main.TStep


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

#######################################################
#PLOTS MEDIAN OF Q_Hii_dot AND SHADES PERCENTILES

median = rag.median_y_values(len(A_P_L_Lya),main.Q_Hii_dot(z, A_P_L_Lya, B_P_L_Lya, C_P_L_Lya, A_f_esc, B_f_esc,Q_Hii))
plt.plot(z,median)

plt.fill_between(z, 0.84*median, 0.16*median)

plt.show()
"""

main.Q_Hii_dot(z,Q,A_P_L_Lya, B_P_L_Lya, C_P_L_Lya, A_f_esc, B_f_esc)
