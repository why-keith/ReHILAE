from scipy.optimize import curve_fit
import numpy as np
import pylab
import matplotlib.pyplot as plt
import math

c_ha=1.36E-12
scale = 13.951808701009465

e1,e2 = 6.67298, 80.96867
f1,f2 = 0.00064, 0.00941
p1,p2,p3,p4,p5,p6 = 1.13869, 39.18853, 0.0013679, -0.0559273, 0.2937793, 26.0952716

e1_err,e2_err = 10.14199, 40.79174
f1_err,f2_err = 0.00013, 0.00364
p1_err,p2_err,p3_err,p4_err,p5_err,p6_err = 0.18238,0.10709,0.0009035,0.0248473,0.2101038,0.5405448

def EW(z, e1=e1, e2=e2):
    return e1*z + e2

def f_esc(z, f1=f1, f2=f2):
    return f1*EW(z) + f2

def rhoLya(z, p1=p1, p2=p2, p3=p3, p4=p4, p5=p5, p6=p6):
    if z <=5.8:
        x = math.log10(1+z)
        return p1*x + p2
    if z>5.8:
        uv = p3*z**3 + p4*z**2 + p5*z + p6
        scaled = uv + scale
        return scaled

def Q_dot(z,P_Lya,fesc,ew):
    return 10**P_Lya / (c_ha*(1-fesc*(0.0042*ew)))

zs = list(np.linspace(0,14,10000))
scipy_fit = []
for z in zs:
    P_Lya = rhoLya(z)#
    fesc = f_esc(z)
    ew = EW(z)
    scipy_fit.append(Q_dot(z,P_Lya,fesc,ew))

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

all_runs = []
for _ in list(range(0,n)):
    all_runs.append([])

for loop in all_runs:
    for z in zs:
        equiWidth_parameter = np.random.randint(2)
        if equiWidth_parameter == 0:
            param = np.random.choice(e1_list)
            ew = EW(z,e1=param)
        else:
            param = np.random.choice(e2_list)
            ew = EW(z,e2=param)
        f_esc_parameter = np.random.randint(2)
        if f_esc_parameter == 0:
            param = np.random.choice(f1_list)
            fesc = f_esc(z,f1=param)
        else:
            param = np.random.choice(f2_list)
            fesc = f_esc(z,f2=param)
        if z <= 5.8:
            P_Lya_parameter = np.random.randint(2)
            if P_Lya_parameter == 0:
                param = np.random.choice(p1_list)
                P_Lya = rhoLya(z,p1=param)
            else:
                param = np.random.choice(p2_list)
                P_Lya = rhoLya(z,p2=param)
        else:
            P_Lya_parameter = np.random.randint(4)
            if P_Lya_parameter == 0:
                param = np.random.choice(p3_list)
                P_Lya = rhoLya(z,p3=param)
            elif P_Lya_parameter == 1:
                param = np.random.choice(p4_list)
                P_Lya = rhoLya(z,p4=param)
            elif P_Lya_parameter == 2:
                param = np.random.choice(p5_list)
                P_Lya = rhoLya(z,p5=param)
            else:
                param = np.random.choice(p6_list)
                P_Lya = rhoLya(z,p6=param)
        loop.append(Q_dot(z, P_Lya, fesc, ew))

median, median_upper_percentile, median_lower_percentile = [],[],[]
for z in zs:
    ind = zs.index(z)
    equiWidths = [i[ind] for i in all_runs]
    median_lower_percentile.append(np.percentile(equiWidths,84))
    median.append(np.median(equiWidths))
    median_upper_percentile.append(np.percentile(equiWidths,16))

plt.figure('Q_dot')
#plt.plot(zs, [math.log10(i) for i in scipy_fit], color='black', label='SciPy fit')
plt.fill_between(zs,  [math.log10(i) for i in median_lower_percentile], [math.log10(i) for i in median_upper_percentile], alpha=0.4, color = "grey", edgecolor = "black", linewidth = 1.2, label=r'68% Confidence Interval')
plt.plot(zs, [math.log10(i) for i in median], "--", label='Median',color='blue')
plt.xlabel('Redshift (z)')
plt.ylabel(r'$\log_{10}(Q_{II} \ [s^{-1} \ Mpc^{-3}])$')
plt.legend()
plt.show()     