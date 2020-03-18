import matplotlib.pyplot as plt
import numpy as np
#from pylab import array
import LAE_Model_new as main
from sys import stdout
import time

start=time.time()

e1,e2 = 6.67298, 80.96867
f1,f2 = 0.00064, 0.00941
p1,p2,p3,p4,p5,p6 = 1.13869, 39.18853, 0.0013679, -0.0559273, 0.2937793, 26.0952716

e1_err,e2_err = 10.14199, 40.79174
f1_err,f2_err = 0.00013, 0.00364
p1_err,p2_err,p3_err,p4_err,p5_err,p6_err = 0.18238,0.10709,0.0009035,0.0248473,0.2101038,0.5405448

n = 100
e1_list = np.random.normal(e1, e1_err*0.2, n)
e2_list = np.random.normal(e2, e2_err*0.2, n)

f1_list = np.random.normal(f1, f1_err*0.2, n)
f2_list = np.random.normal(f2, f2_err*0.2, n)

p1_list = np.random.normal(p1, p1_err*0.2, n)
p2_list = np.random.normal(p2, p2_err*0.2, n)
p3_list = np.random.normal(p3, p3_err*0.2, n)
p4_list = np.random.normal(p4, p4_err*0.2, n)
p5_list = np.random.normal(p5, p5_err*0.2, n)
p6_list = np.random.normal(p6, p6_err*0.2, n)

all_runs = []
ts = np.linspace(0.051,14,100000) # time in Gyr
zs= list(((((28./(ts))-1.)**(1./2.)-1.))) # conversion from Gyr to redshift

for _ in list(range(0,n)):
    equiWidth_parameter = np.random.randint(2)
    if equiWidth_parameter == 0:
        param = np.random.choice(e1_list)
        ew = [param, e2]
    else:
        param = np.random.choice(e2_list)
        ew = [e1, param]
    f_esc_parameter = np.random.randint(2)
    if f_esc_parameter == 0:
        param = np.random.choice(f1_list)
        fesc = [param, f2]
    else:
        param = np.random.choice(f2_list)
        fesc = [f1, param]
    P_Lya_parameter_power = np.random.randint(2)
    if P_Lya_parameter_power == 0:
        param = np.random.choice(p1_list)
        P_Lya = [param,p2,p3,p4,p5,p6]
    else:
        param = np.random.choice(p2_list)
        P_Lya = [p1,param,p3,p4,p5,p6]
    P_Lya_parameter_cubic = np.random.randint(4)
    if P_Lya_parameter_cubic == 0:
        param = np.random.choice(p3_list)
        P_Lya[2] = param
    elif P_Lya_parameter_cubic == 1:
        param = np.random.choice(p4_list)
        P_Lya[3] = param
    elif P_Lya_parameter_cubic == 2:
        param = np.random.choice(p5_list)
        P_Lya[4] = param
    else:
        param = np.random.choice(p6_list)
        P_Lya[5] = param
    parameters = (ew,fesc,P_Lya)
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