import matplotlib.pyplot as plt
import numpy as np
#from pylab import array
import LAE_Model_new as main
######################################################## PARAMETERS
e1,e2 = 6.67298, 80.96867
f1,f2 = 0.00064, 0.00941
p1,p2,p3,p4,p5,p6 = 1.13869, 39.18853, 0.0013679, -0.0559273, 0.2937793, 26.0952716

e1_err,e2_err = 10.14199, 40.79174
f1_err,f2_err = 0.00013, 0.00364
p1_err,p2_err,p3_err,p4_err,p5_err,p6_err = 0.18238,0.10709,0.0009035,0.0248473,0.2101038,0.5405448

######################################################## SAMPLE RANDOM PARAMETERS
n = 1000
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

######################################################## SIMULATION
all_runs = []
ts = np.linspace(0.051,14,100000) 
zs= list(((((28./(ts))-1.)**(1./2.)-1.)))

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
    fraction = [i[ind] for i in all_runs]
    median_lower_percentile.append(np.percentile(fraction,84))
    median.append(np.median(fraction))
    median_upper_percentile.append(np.percentile(fraction,16))

######################################################## DATA POINTS
mason_x, mason_y = [6.996610169491525], [1-0.5910931174089071]
mason_xerr = [[0.50169491525], [0.50169491525]]
mason_yerr = [[0.15060728744], [0.10931174089]]
lya_fraction_x, lya_fraction_y = [6.894915254237289], [1-0.47044534412955463]
LAE_clustering_x, LAE_clustering_y = [6.596610169491526], [1-0.42995951417004075]
LAE_err = [[0.19661016949], [3.6]]
dark_fraction_x, dark_fraction_y = [6.094915254237288,5.898305084745763,5.593220338983051], [1-0.322267206477733,1-0.04210526315789509,1-0.022672064777328416]
QSO_damping_wings_x, QSO_damping_wings_y = [7.532203389830496,7.098305084745764], [1-0.5497975708502038,1-0.40000000000000013]
QSO_yerr = [[0.1991902834,0.1991902834],[0.2,0.2]]
Planck_x, Planck_y = [8.298305084745763], [1-0.5004048582995955]

######################################################## PLOTS
plt.figure('All iterations')
for line in all_runs:
    plt.plot(zs, line)
plt.xlabel('Redshift (z)')
plt.ylabel(r'Fraction of Ionised Hydrogen ($Q_{II}$)')
plt.axvspan(6, 10, color = "lightgrey", alpha = 0.3, edgecolor = "black", linewidth = 5)

plt.figure('Fraction of Ionised Hydrogen LAE')
plt.fill_between(zs,median_lower_percentile,median_upper_percentile,alpha=0.4,color ="steelblue",edgecolor="black",linewidth = 1.2,zorder=3)
plt.plot(zs, median, color='black',zorder=3)
plt.axvspan(6,10,color ="lightgrey",alpha=0.4,edgecolor="black",linewidth=5,zorder=1)
plt.scatter(mason_x,mason_y,edgecolors='red',facecolors='red',marker='*',s=250,label='Mason+2018',zorder=10)
plt.errorbar(mason_x,mason_y,yerr=mason_yerr,xerr=mason_xerr,capsize=5,color='red',zorder=10)
plt.scatter(lya_fraction_x,lya_fraction_y,edgecolors='black',facecolors='white',marker='*',s=250,label=r'Ly$\alpha$ fraction',zorder=10)
plt.scatter(LAE_clustering_x,LAE_clustering_y,facecolors='black',marker='s',s=50,label='LAE clustering',zorder=10)
plt.errorbar(LAE_clustering_x,LAE_clustering_y,xerr=LAE_err,capsize=5,color='black',zorder=10)
plt.scatter(dark_fraction_x,dark_fraction_y,facecolors='black',marker='o',s=50,label='Dark fraction',zorder=10)
plt.scatter(QSO_damping_wings_x,QSO_damping_wings_y,facecolors='black',marker='D',s=50,label='QSO daming wings',zorder=10)
plt.errorbar(QSO_damping_wings_x,QSO_damping_wings_y,yerr=QSO_yerr,capsize=5,color='black',ls='none',zorder=10)
plt.scatter(Planck_x,Planck_y,edgecolors='black',facecolors='black',marker='p',s=50,label='Planck 2016',zorder=10)
plt.xlabel('Redshift (z)')
plt.ylabel(r'Fraction of Ionised Hydrogen ($Q_{II}$)')
plt.legend()

plt.show()