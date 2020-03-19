import matplotlib.pyplot as plt
import numpy as np

data=np.load("UV_saves/UV_C=1_xi_constant.npy")
zs, median, median_lower_percentile, median_upper_percentile=data[0],data[1],data[2],data[3]

data3=np.load("UV_saves/UV_C=3_xi_constant.npy")
zs3, median3, median_lower_percentile3, median_upper_percentile3=data3[0],data3[1],data3[2],data3[3]

data10=np.load("UV_saves/UV_C=10_xi_constant.npy")
zs10, median10, median_lower_percentile10, median_upper_percentile10=data10[0],data10[1],data10[2],data10[3]

plt.figure('Fraction of Ionised Hydrogen LAE')
plt.fill_between(zs, median_lower_percentile,  median_upper_percentile, alpha=0.4,  edgecolor = "black", linewidth = 1.2, label=r'68% Confidence Interval, C=1')
plt.fill_between(zs3, median_lower_percentile3,  median_upper_percentile3, alpha=0.4,  edgecolor = "blue", linewidth = 1.2, label=r'68% Confidence Interval, C=3')
plt.fill_between(zs10, median_lower_percentile10,  median_upper_percentile10, alpha=0.4,  edgecolor = "red", linewidth = 1.2, label=r'68% Confidence Interval, C=10')
plt.plot(zs, median, label='C=1')
plt.plot(zs, median, label='C=3')
plt.plot(zs, median, label='C=10',color="black")
plt.axvspan(6, 10, color = "lightgrey", alpha = 0.4, edgecolor = "black", linewidth = 5)
plt.xlabel('Redshift (z)')
plt.ylabel(r'Fraction of Ionised Hydrogen ($Q_{H_{II}}$)')
plt.legend()

plt.show()