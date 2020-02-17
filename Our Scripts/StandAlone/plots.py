import math
import numpy as np
import pylab
import matplotlib.pyplot as plt
import LyCmodel as model 

def EW_z():
    """
    Data from "Table_C3_Calhau19_Stacking_LAEs_X_rays_v1.fits"
    """
    redshift = [2.5, 2.8, 2.9, 3.1, 3.3, 3.7, 4.1, 4.5, 4.8, 5.0, 5.3]
    equiWidth = [117.134349876, 122.113769531, 172.679885864, 143.5936203, 149.371124268, 96.3648490905, 189.002449036, 181.127731323, 95.0448436864, 125.935153962, 143.343048096]
    x = pylab.array(redshift)
    y = pylab.array(equiWidth)
    p2 = pylab.polyfit(x, y, 1.0)
    p = pylab.poly1d(p2)

    ts = np.linspace(0.051,14,100000) # time in Gyr
    zs= ((((28./(ts))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift
    EW = [p(i) for i in zs]

    plt.plot(zs,EW,color='red')
    plt.plot(redshift, equiWidth,'o',color='black',markersize=4)
    plt.xlabel('Redshift ($z$)')
    plt.ylabel(r'Equivalent Width ($EW_{0}$) [$\mathring{A}$]')
    plt.title('Plot to Show $EW_{0}$ Against Redshift')
    plt.show()

def EW_Fesc():
    """
    Data from https://www.aanda.org/articles/aa/pdf/2017/01/aa29264-16.pdf
    """
    EW = [79, 129, 83, 98, 75, 29, 15, 4]
    escape = [0.132, 0.074, 0.072, 0.058, 0.056, 0.045, 0.032, 0.01]

    x = pylab.array(EW)
    y = pylab.array(escape)
    p1 = pylab.polyfit(x,y,1.0)
    p = pylab.poly1d(p1)
    z = np.linspace(0,150,1000)
    f = [p(i) for i in z]

    EW_ave = [140.321828866]
    LyC_f = [p(140.321828866)]

    plt.plot(z,f,color='red')
    plt.plot(EW, escape,'o',color='black',markersize=4)
    plt.plot(EW_ave, LyC_f, 'o',color='blue')
    plt.xlabel(r'Equivalent Width ($EW_{0}$) [$\mathring{A}$]')
    plt.ylabel('LyC Escape Fraction')
    plt.title('Plot to Show $f_{esc, LyC}$ against $EW_{0}$')
    plt.show()

def PLya():
    """
    Data from Table_C2_Sobral18_SSC4K_compilation.fits
    """
    redshift = [2.2,2.5,2.8,3.0,3.2, 3.3, 3.7, 4.1, 4.6, 4.8, 5.1, 5.3, 5.8 ]
    P_Lya = [0.52, 0.74, 0.77, 0.88, 0.84, 0.85, 1.01, 0.87, 1.19, 1.12, 1.27, 1.08, 1.10]
    x = pylab.array(redshift)
    y = pylab.array(P_Lya)  # data from SC4K Sobral 
    p2 = pylab.polyfit(x, y, 2.0)
    p = pylab.poly1d(p2)

    ts = np.linspace(0.051,14,100000) # time in Gyr
    zs= ((((28./(ts))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift
    rho = []
    for i in zs:
        if p(i)<0:
            rho.append(0)
            continue
        #peak=np.max(y)
        #peak_position = np.where(y==peak)
        #cutoff=x[peak_position[0][0]]
        #if i > cutoff:
        #    rho.append(p(cutoff) * 10**40)
        else:
            rho.append(p(i)*10**40)
    L_Lya = [i*10**(40) for i in P_Lya]

    plt.plot(zs,rho,color='red')
    plt.plot(redshift, L_Lya,'o',color='black',markersize=4)
    plt.xlabel('Redshift ($z$)')
    plt.ylabel(r'Ly$\alpha$ Luminosity ($L_{Ly_{\alpha}}$)')
    plt.title(r'Plot to Show the Ly$\alpha$ Luminosity Against Redshift')
    plt.show()

def n_ion():
    ts = np.linspace(0.051,14,100000) # time in Gyr
    zs= ((((28./(ts))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift

    n_ion = []
    for i in zs:
        n_ion.append(model.n_ion_dot_LyC(i))

    plt.plot(zs, n_ion)
    plt.xlabel('Redshift ($z$)')
    plt.ylabel(r'$\dot{n}_{ion} \ [cm^{-3} \ s^{-1}]$')
    plt.title(r'Plot to Show $\dot{n}_{ion}$ Against Redshift')
    plt.show()

def Qion_Lyc():
    ts = np.linspace(0.051,14,10000000) # time in Gyr
    zs= ((((28./(ts))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift

    q_ii = []
    for i in zs:
        q_ii.append(model.Q_ion_LyC(i))

    plt.plot(zs, q_ii)
    plt.xlabel('Redshift ($z$)')
    plt.ylabel(r'$\dot{Q}_{ion, Ly\alpha} \ [Mpc^{-3} \ s^{-1}]$')
    plt.title(r'Plot to Show $\dot{Q}_{ion, Ly\alpha}$ Against Redshift')
    plt.show()

def t_rec():
    ts = np.linspace(0.051,14,100000) # time in Gyr
    zs= ((((28./(ts))-1.)**(1./2.)-1.)) # conversion from Gyr to redshift

    recTime = []
    for i in zs:
        recTime.append(model.t_rec(i))
    
    plt.plot(zs, recTime)
    plt.xlabel('Redshift ($z$)')
    plt.ylabel(r'Recombination Time ($t_{rec}) \ [s]$')
    plt.title(r'Plot to Show $t_{rec}$ Against Redshift')
    plt.show()

PLya()