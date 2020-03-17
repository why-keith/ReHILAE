import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes,InsetPosition,mark_inset)

"""
see https://scipython.com/blog/inset-plots-in-matplotlib/ for reference
"""

x = np.linspace(0,2*np.pi,10000) # generate bottom x-axis
y = [np.sin(angle) for angle in x] # y values that both x-axes will share
x_prime = [np.degrees(angle) for angle in x] # generate top y-axis

fig, ax = plt.subplots()
ax.plot(x,y,label=r'$\sin(\theta)$') # plot bottom x and y data
ax.set_xlabel(r'$\theta$ (radians)')
ax1 = plt.twiny() # create instance for top axis
ax1.set_xlabel(r'$\theta$ (degrees)')
ax1.plot(x_prime,y) # plot top x and y data
ax2 = plt.axes([0,0,1,1]) # create inset instance, no idea what the values do
ip = InsetPosition(ax1, [0.2,0.1,0.3,0.3]) # all of these are fractional, format: x coordinate of plot. y coordinate of plot, height, width
ax2.set_axes_locator(ip) # assign inset position to inset instance
mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5') # draw grey lines from inset position to data on plot
ax2_x, ax2_y = [],[]
for (x,y) in zip(x_prime,y):
    if 89.5<x<90.5: # the data we are focusing on
        ax2_x.append(x)
        ax2_y.append(y)
ax2.plot(ax2_x,ax2_y) # plot the data on the inset
ax.legend()
plt.show()