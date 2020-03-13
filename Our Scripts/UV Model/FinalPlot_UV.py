import matplotlib.pyplot as plt
import numpy as np
import math
import random_array_generator as rag
import pandas as pd
import Main_UV

startT = Main_UV.startT
finishT = Main_UV.finishT
TStep = Main_UV.TStep

P1_1 = [-0.0184311]
P1_1_error = [0.0023189]

P1_2 = [-0.0185036]
P1_2_error = [0.0460448]

P1_3 = [26.8715339]
P1_3_error = [0.1969757]

P1_1_P_UV = rag.random_Arrays(len(P1_1),P1_1,P1_1_error,P1_1_error)
P1_2_P_UV = rag.random_Arrays(len(P1_2),P1_2 ,P1_2_error,P1_2_error)
P1_3_P_UV = rag.random_Arrays(len(P1_3),P1_3 ,P1_3_error,P1_3_error)

data = []
z,t = Main_UV.redshift(startT, finishT, TStep)
for i,j,k in zip(P1_1_P_UV,P1_2_P_UV,P1_3_P_UV):
    arguements = (i[0],j[0],k[0])
    data.append((Main_UV.main_UV(arguements))) # TODO write this return to file and then plot after loop




#######################################################
#PLOTS MEDIAN OF Q_Hii_dot AND SHADES PERCENTILES
median, median_lower_percentile, median_upper_percentile = rag.median_y_values(len(data[0]),data)
plt.figure('Ionised_Hydrogen_UV  10000')
plt.title("Fraction of ionised H with respect to redshift")
plt.xlabel("Redshift (z)")
plt.ylabel("Fractions of Ionised Hydrogen")
plt.plot(z,median, color = "steelblue")
#plt.legend()

plt.fill_between(z,  median_lower_percentile, median_upper_percentile, alpha=0.4, color = "lightblue", edgecolor = "black", linewidth = 1.2)
plt.fill_betweenx(median,6,10, color = "lightgrey", alpha = 0.3, edgecolor = "black", linewidth = 1.2)
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
plt.show()          
    
    
df=pd.DataFrame(data=Main_UV.param_dict)
try:
    pd.DataFrame.to_csv(df,"Inital_parameters.csv")
except:
    print("Writing to\"Inital_parameters.csv\" failed")