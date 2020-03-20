import matplotlib.pyplot as plt
import matplotlib
import numpy as np

quad = np.load('LAE_results_quad.npz')
quad_zs = quad['redshift']
quad_lower = quad['lower']
quad_median = quad['median']
quad_higher = quad['upper']

hybrid = np.load('LAE_results.npz')
hybrid_zs = hybrid['redshift']
hybrid_lower = hybrid['lower']
hybrid_median = hybrid['median']
hybrid_higher = hybrid['upper']

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

plt.figure('LAE Models Comparison')

plt.fill_between(hybrid_zs,hybrid_lower,hybrid_higher,alpha=0.4,color ="steelblue",edgecolor="black",linewidth = 1.2,label=r'$1\sigma$ conf. interval (hybrid)',zorder=3)
plt.plot(hybrid_zs, hybrid_median, color='black',label=r'ReHiLAE (hybrid $\rho_{Ly\alpha}$)',zorder=3)

plt.fill_between(quad_zs,quad_lower,quad_higher,alpha=0.4,color ="green",edgecolor="black",linewidth = 1.2,label=r'$1\sigma$ conf. interval (quadratic)',zorder=3)
plt.plot(quad_zs, quad_median, linestyle='--',color='black',label=r'ReHiLAE (quadratic $\rho_{Ly\alpha}$)',zorder=3)


plt.axvspan(6,10,facecolor ="lightgrey",alpha=0.4,edgecolor="grey",linewidth=5,zorder=1)

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
plt.xlabel('Redshift (z)')
plt.ylabel(r'Fraction of Ionised Hydrogen ($Q_{HII}$)')
plt.legend()

matplotlib.rcParams['lines.linewidth'] = 6
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