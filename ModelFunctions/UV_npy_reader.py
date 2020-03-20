import matplotlib.pyplot as plt
import matplotlib
import numpy as np

data=np.load("UV_saves/UV_C=1_xi_constant.npy")
zs, median, median_lower_percentile, median_upper_percentile=data[0],data[1],data[2],data[3]

#data3=np.load("UV_saves/UV_C=3_xi_variable.npy")
#zs3, median3, median_lower_percentile3, median_upper_percentile3=data3[0],data3[1],data3[2],data3[3]

#data10=np.load("UV_saves/UV_C=10_xi_constant.npy")
#zs10, median10, median_lower_percentile10, median_upper_percentile10=data10[0],data10[1],data10[2],data10[3]


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


plt.figure('Fraction of Ionised Hydrogen LAE')
#plt.fill_between(zs, median_lower_percentile,  median_upper_percentile, alpha=0.4,  edgecolor = "black", linewidth = 1.2, label=r'$1\sigma$ conf. interval, C=1')
plt.fill_between(zs,median_lower_percentile,median_upper_percentile,alpha=0.4,color ="steelblue",edgecolor="black",linewidth = 1.2,label=r'$1\sigma$ conf. interval',zorder=3)
#plt.fill_between(zs10, median_lower_percentile10,  median_upper_percentile10, alpha=0.4,  edgecolor = "red", linewidth = 1.2, label=r'$1\sigma$ conf. interval, C=10')
#plt.plot(zs, median, label='C=1')
plt.plot(zs, median, color='black',label='median',zorder=3)
#plt.plot(zs, median10, label='C=10',color="black")
plt.axvspan(6, 10, color = "lightgrey", alpha = 0.4, edgecolor = "black", linewidth = 5)
plt.xlabel('Redshift (z)')
plt.ylabel(r'Fraction of Ionised Hydrogen ($Q_{H_{II}}$)')

plt.scatter(mason_x,mason_y,edgecolors='red',facecolors='red',marker='*',s=250,label='Mason+2018',zorder=10)
plt.errorbar(mason_x,mason_y,yerr=mason_yerr,xerr=mason_xerr,capsize=5,color='red',zorder=10)
plt.scatter(lya_fraction_x,lya_fraction_y,edgecolors='black',facecolors='white',marker='*',s=250,label=r'Ly$\alpha$ fraction',zorder=10)
plt.scatter(LAE_clustering_x,LAE_clustering_y,facecolors='black',marker='s',s=50,label='LAE clustering',zorder=10)
plt.errorbar(LAE_clustering_x,LAE_clustering_y,xerr=LAE_err,capsize=5,color='black',zorder=10)
plt.scatter(dark_fraction_x,dark_fraction_y,facecolors='black',marker='o',s=50,label='Dark fraction',zorder=10)
plt.scatter(QSO_damping_wings_x,QSO_damping_wings_y,facecolors='black',marker='D',s=50,label='QSO damping wings',zorder=10)
plt.errorbar(QSO_damping_wings_x,QSO_damping_wings_y,yerr=QSO_yerr,capsize=5,color='black',ls='none',zorder=10)
plt.scatter(Planck_x,Planck_y,edgecolors='black',facecolors='black',marker='p',s=50,label='Planck 2016',zorder=10)

plt.tick_params(which='both',direction='in',right=True,top=True)

plt.xlim(0,20)

plt.legend()

matplotlib.rcParams['lines.linewidth'] = 1.2
matplotlib.rcParams['axes.linewidth'] = 2.0
matplotlib.rcParams['xtick.major.size'] = 9
matplotlib.rcParams['xtick.minor.size'] = 5
matplotlib.rcParams['xtick.major.width'] = 1.9
matplotlib.rcParams['xtick.minor.width'] = 1.3
matplotlib.rcParams['ytick.major.size'] = 9
matplotlib.rcParams['ytick.minor.size'] = 4
matplotlib.rcParams['ytick.major.width'] = 1.9
matplotlib.rcParams['ytick.minor.width'] = 1.3

plt.show()