import matplotlib.pyplot as plt
import numpy as np
#from pylab import array
import UV_Model_new as main
from sys import stdout
import time

start=time.time()

C1 = 0.0013679194302549992
C1_err = 0.0009034769902604823

C2 = -0.05592731916980156
C2_err = 0.024847250273459104

C3 = 0.29377928649363305
C3_err = 0.21010382042922288

C4 = 26.095271603086044
C4_err = 0.5405447975405038

f1=2.3/100
f1_err=0.1/100
f2=1.17
f2_err=0.02


n = 1000


C1_list = np.random.normal(C1, C1_err*0.2, n)
C2_list = np.random.normal(C2, C2_err*0.2, n)
C3_list = np.random.normal(C3, C3_err*0.2, n)
C4_list = np.random.normal(C4, C4_err*0.2, n)

f1_list = np.random.normal(f1, f1_err*0.2, n)
f2_list = np.random.normal(f2, f2_err*0.2, n)

"""
C1_list = np.random.normal(C1, C1_err*0.2, n)
c2_list = np.random.normal(c2, c2_err*0.2, n)

C3_list = np.random.normal(C3, C3_err*0.2, n)
C4_list = np.random.normal(C4, C4_err*0.2, n)

p1_list = np.random.normal(p1, p1_err*0.2, n)
p2_list = np.random.normal(p2, p2_err*0.2, n)
p3_list = np.random.normal(p3, p3_err*0.2, n)
p4_list = np.random.normal(p4, p4_err*0.2, n)
p5_list = np.random.normal(p5, p5_err*0.2, n)
p6_list = np.random.normal(p6, p6_err*0.2, n)
"""
all_runs = []
ts = np.linspace(0.051,14,100000) # time in Gyr
zs= list(((((28./(ts))-1.)**(1./2.)-1.))) # conversion from Gyr to redshift

for _ in list(range(0,n)):   
    
    percent=100*(_/n)   
    stdout.write("\rRandomising parameters - {}%       ".format(round(percent,3)))
    stdout.flush()
    
    
    
    c_parameter=np.random.randint(4)
    if c_parameter==0:
        param=np.random.choice(C1_list)
        c=[param,C2,C3,C4]
    elif c_parameter==1:
        param=np.random.choice(C2_list)
        c=[C1,param,C3,C4]
        
    elif c_parameter==2:
        param=np.random.choice(C3_list)
        c=[C1,C2,param,C4]
        
    else:
        param=np.random.choice(C4_list)
        c=[C1,C2,C3,param]
        
        
    f_parameter=np.random.randint(2)
    if f_parameter==0:
        param=np.random.choice(f1_list)
        f=[param,f2]        
    else:
        param=np.random.choice(f2_list)
        f=[f1,param]
        

        
    parameters = (c,f)
    arguements = [i for sub in parameters for i in sub]
    all_runs.append(main.main(ts,arguements))

median, median_upper_percentile, median_lower_percentile = [],[],[]
for z in zs:
    ind = zs.index(z)
    
    percent=100*ind/len(zs)    
    stdout.write("\rGenerating median arrays - {}%       ".format(round(percent,3)))
    stdout.flush()
  
    
    
    fraction = [i[ind] for i in all_runs]
    median_lower_percentile.append(np.percentile(fraction,84))
    median.append(np.median(fraction))
    median_upper_percentile.append(np.percentile(fraction,16))
    
Data=np.array([zs, median, median_lower_percentile, median_upper_percentile])
path="UV_saves/UV_C="+str(main.C)+"_xi_variable"
np.save(path,Data)
print("\nData saved")

plt.figure('All iterations')
for line in all_runs:
    plt.plot(zs, line)
plt.xlabel('Redshift (z)')
plt.ylabel(r'Fraction of Ionised Hydrogen ($Q_{II}$)')
plt.axvspan(6, 10, color = "lightgrey", alpha = 0.3, edgecolor = "black", linewidth = 5)

plt.figure('Fraction of Ionised Hydrogen LAE')
plt.fill_between(zs, median_lower_percentile,  median_upper_percentile, alpha=0.4, color = "steelblue", edgecolor = "black", linewidth = 1.2, label=r'68% Confidence Interval')
plt.plot(zs, median, label='Median',color='black')
plt.axvspan(6, 10, color = "lightgrey", alpha = 0.4, edgecolor = "black", linewidth = 5)
plt.xlabel('Redshift (z)')
plt.ylabel(r'Fraction of Ionised Hydrogen ($Q_{II}$)')
plt.legend()

plt.show()



print("\nTime elapsed = {}s".format(round(time.time()-start,2)))