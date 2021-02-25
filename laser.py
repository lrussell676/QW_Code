# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:42:39 2021

@author: Lewis
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [12, 8]
plt.rcParams['figure.dpi'] = 100
plt.rcParams.update({'font.size': 16})
import scipy.optimize as sc

"""

04/02/2021, 11/02/2021, 13/02/2021

"""
linspace_size = 100
h_2 = (6.62607015e-34)**2
h_bar2 = ((6.62607015e-34)/(2*np.pi))**2
m = 0.06*(9.11e-31)
a = np.linspace(0.1e-9,10e-9,linspace_size)
al_ratio = np.linspace(0, 0.4, 100000)
E_0 = ((h_2/(8*m*(a**2))))/(1.6e-19)
E_g = (1.424 + 1.247*al_ratio)

Figure, Plot4 = plt.subplots(1,1)
Figure.patch.set_facecolor('xkcd:mint green')
Figure.suptitle("Al Ratio Plot")
Plot4.set(ylabel = "$E_{g}$ (eV)", xlabel = "AlGa Fraction (x)")
Plot4.plot(al_ratio,E_g)
Plot4.plot(al_ratio[83514],E_g[83514], marker='x', markersize=15)
Plot4.grid(b=True, which='major', linestyle=':', linewidth='2')
Plot4.legend(['$E_{g}$','X Marker Point of E_pump','X Marker Point of E_emission'],loc='best')
plt.show()

h = (6.62607015e-34)
E_pump = (h* ((3e8)/(675e-9)) )/(1.6e-19)
E_emission = (h* ((3e8)/(815e-9)) )/(1.6e-19)

print(f"\n \n E_pump is {E_g[83514]}eV at AlGa Fraction {al_ratio[83514]}")

print("\n Take 'A' Barrier Composition as Al_{0.330}Ga_{0.670}As")

print("\n Take 'B' Well Composition as GaAs")

Figure, Plot5 = plt.subplots(1,1)
Figure.patch.set_facecolor('xkcd:mint green')
Figure.suptitle("$E_{0}$ vs Well Width")
Plot5.set(ylabel = "$E_{0}$/eV", xlabel = "a (m)")
Plot5.plot(a,E_0)
Plot5.grid(b=True, which='major', linestyle=':', linewidth='2')


m = (9.11e-31)
me = 0.06*m
mhh = 0.57*m
a = np.linspace(0.1e-9,7e-9,linspace_size)

E_0e = ((h_2/(8*me*(a**2))))/(1.6e-19)
E_0hh = ((h_2/(8*mhh*(a**2))))/(1.6e-19)

Figure, Plot6 = plt.subplots(1,1)
Figure.suptitle("$E_{0}$ vs Well Width")
Plot6.set(ylabel = "$E_{0}$/eV", xlabel = "a (nm)")
Plot6.plot(a*(1e9),E_0e, label = '$E_{e(0)}$ =' r'$\frac{h^2}{8m_{0}a^2}$')
Plot6.plot(a*(1e9),E_0hh, label = '$E_{hh(0)}$ =' r'$\frac{h^2}{8m_{hh}a^2}$')
Plot6.legend(loc = 'best')
Plot6.grid(b=True, which='major', linestyle=':', linewidth='2')
plt.minorticks_on()
Plot6.grid(b=True, which='minor', linestyle='-', linewidth='0.2')

plt.ylim(0, 5)
plt.show()

h = 6.62607004e-34
c = 299792458
m0 = 9.11e-31

V0 = (1.247*0.33)*(1.66e-19)
V0e = 2/3*V0
V0h = 1/3*V0
a = np.arange(1e-9,11e-9,0.1e-9)

Figure, (SubPlot1,SubPlot2) = plt.subplots(1,2,constrained_layout=True)

E0e_main, E0h_main, zi_e, zi_h = np.empty(len(a)), \
                    np.empty(len(a)), \
                    np.empty(len(a)), \
                    np.empty(len(a))
c = 0
for i in a:
    E0e = (h**2)/(8*me*(i**2))
    E0h = (h**2)/(8*mhh*(i**2))
    
    E0e_main[c] = E0e
    E0h_main[c] = E0h

    de = V0e/E0e
    dh = V0h/E0h
    
    zi = np.linspace(0,20,1000)
    ye = np.sqrt(((de*(np.pi**2))/4) - zi**2)
    yh = np.sqrt(((dh*(np.pi**2))/4) - zi**2)
    y2 = zi*(np.tan(zi))

    ye_new = ye.astype(np.float)
    yh_new = yh.astype(np.float)
    y2new = y2.astype(np.float)

    ye_new[ye <= 0] = np.nan
    yh_new[yh <= 0] = np.nan
    y2new[y2 <= 0] = np.nan


    SubPlot1.plot(zi,ye_new)
    SubPlot1.plot(zi,y2new)
    SubPlot2.plot(zi,yh_new)
    SubPlot2.plot(zi,y2new)
    
    if c > 0 and c <= 10:
        EstIntElectron = np.array([0.3])
    elif c > 10 and c <= 20:
        EstIntElectron = np.array([0.5])
    elif c > 20 and c <=30:
        EstIntElectron = np.array([0.8])
    else:
        EstIntElectron = np.array([1.2])
    
    if c > 0 and c <= 10:
        EstIntHole = np.array([0.5])
    elif c > 10 and c <= 20:
        EstIntHole = np.array([1])
    else:
        EstIntHole = np.array([1.55])
        
    def funcelectron(x):
        return np.sqrt(((de*(np.pi**2))/4) - x**2) - x*(np.tan(x))

    def funchole(x):
        return np.sqrt(((dh*(np.pi**2))/4) - x**2) - x*(np.tan(x))

    root_e = sc.fsolve(funcelectron,EstIntElectron)
    root_h = sc.fsolve(funchole,EstIntHole)

    zi_e[c] = root_e
    zi_h[c] = root_h
    c += 1
      
E1e = ((4*zi_e**2)/(np.pi**2))*E0e_main
E1h = ((4*zi_h**2)/(np.pi**2))*E0h_main
    
SubPlot1.set(xlabel = '$\u03BE$', \
             ylabel = 'q($\u03BE$), p($\u03BE$)', \
             title = 'q(\u03BE),p(\u03BE) vs \u03BE For effective electron mass', \
             ylim = ([0,5]), xlim = ([0,2]))
SubPlot2.set(xlabel = '$\u03BE$', \
             ylabel = 'q($\u03BE$), p($\u03BE$)', \
             title = 'q(\u03BE),p(\u03BE) vs \u03BE For effective hole mass', \
             ylim = ([0,8]), xlim = ([0,2]))
SubPlot1.legend(['q(\u03BE) = $\sqrt{d_{electron}\pi^{2}/4 - \u03BE^2}$', \
                 'p(\u03BE) = $\u03BE$tan($\u03BE$)'], loc = 'best')
SubPlot1.grid(True)
SubPlot2.grid(True)
plt.show()

#Plotting Transition Energy vs a
a = a/1e-9
Et = (E1e + E1h + 1.424*1.66e-19)/(1.66e-19)
Figure, SubPlot1 = plt.subplots(1,1,constrained_layout=True)
SubPlot1.plot(a,Et)
SubPlot1.axhline(y=E_emission, linewidth='0.8', \
                 linestyle='--', color='orange')
SubPlot1.axvline(x=5.22, linewidth='0.8', linestyle='--', \
                 color='orange')
SubPlot1.set(xlabel = 'a (nm)', \
             ylabel = '$E_{t}$ (eV)',\
                 title = '$E_{t}$ vs Well Width')
SubPlot1.legend([r'$E_{t}$ = $E_{1e} + E_{1h} + V_{0}$',\
                 '$E_{emission}$ at a=5.22nm'], \
                loc = 'best', prop={'size': 15})
SubPlot1.grid(b=True, which='major', linestyle=':', linewidth='2')
plt.minorticks_on()
SubPlot1.grid(b=True, which='minor', linestyle='-', linewidth='0.2')
plt.xlim(3,8)
plt.ylim(1.475,1.575)
plt.show()
    