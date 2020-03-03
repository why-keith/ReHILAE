import numpy as np
import matplotlib
import math
from astropy.io import fits 
import matplotlib.pyplot as plt
import copy

iterations = 100


############################################
# Return condicence levels
def plot_confidence_Lya(x_in,ysample,colour='g'):
    lower = np.percentile(ysample,16,axis=0)  # Equivalent to 1 sigma
    upper = np.percentile(ysample,84,axis=0) # Equivalent to 1 sigma
    lower2 = np.percentile(ysample,100-97.72,axis=0)  #, Equivalent to 2 sigma
    upper2 = np.percentile(ysample,97.72,axis=0) # Equivalent to 2 sigma
    lower3 = np.percentile(ysample,0.15,axis=0)  #, Equivalent to 3 sigma
    upper3 = np.percentile(ysample,99.85,axis=0) # Equivalent to 3 sigma
    return x_in,lower,upper, lower2, upper2, lower3, upper3

####################################
# Gaussian Function
def gaussian(x,x0,sigma):
    return np.exp(-(x-x0)**2/(2*sigma**2))

####################################
# Double gaussian, asymmetric
# phi_dummy generates the space of values to pick from
# phi_probs generates the PDF for a 1sigma error around each phi
def double_normal(phi,err_down,err_up,size):
   phi_dummy = np.linspace(phi-5.*err_down,phi+5.*err_up,100000) #generate list from -5 to +5 error
   phi_probs = np.append(gaussian(phi_dummy[phi_dummy<phi],phi,err_down),gaussian(phi_dummy[phi_dummy>=phi],phi,err_up))
   return np.random.choice(phi_dummy,size,p=phi_probs/np.sum(phi_probs))
   
##################################
# Generates random error arrays

def random_Arrays(x,y,error_down_y,error_up_y):
    master_list = []
    
    for i in range(iterations): # Creates 'iteration' number of arrays with new random guassian distributed data points
        y_new_list = [0 for i in range(len(y))]
        for j in range(len(y)): # Generates an array with random guassian distributed data points
            y_new_list[j]=(double_normal(y[j],error_down_y[j],error_up_y[j],1)[0])
                
        master_list.append(y_new_list) # Appends new array to master_list
    return master_list
    
#################################
# Generates median y values

def median_y_values(length_of_each_array,array_of_random_arrays):
    median_y_array = []
    upper_percentile = []
    lower_percentile = []
    for i in range(length_of_each_array):
        Y = []
        for j in range(len(array_of_random_arrays)):
            Y.append(array_of_random_arrays[j][i])
        median_y_array.append(np.median(Y))
        upper_percentile.append(np.percentile(Y,84))
        lower_percentile.append(np.percentile(Y,16))
        
    #print(median_y_array)

    return median_y_array, upper_percentile, lower_percentile
    
    #plt.plot(x,median_y_array)
    #plt.show()
"""    
##################################
# Finds the max and min of the generated arrays
    max_list = copy.deepcopy(y)
    min_list = copy.deepcopy(y)
    
    for i in range(len(master_list[0])): # Determines the max and min of the random arrays in master_list
        for j in range(len(master_list)):
            if master_list[j][i] > max_list[i]:
                max_list[i] = master_list[j][i]
            if master_list[j][i] < min_list[i]:
                min_list[i] = master_list[j][i]
    return min_list, max_list
################################
# Plots a graph with the max and min error arrays and fills inbetween
   
plt.fill_between(x,min_list,max_list,lw=1,color='#0066ff',alpha=0.1,zorder = 90)
    
plt.plot(x,y)
plt.show()
#################################
"""