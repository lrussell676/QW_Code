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

18/02/2021

"""

linspace_size = 300
c = (2.998e8)
eV = (1.602e-19)
h = 6.626e-34
h_2 = (6.626e-34)**2
h_bar2 = ((6.626e-34)/(2*np.pi))**2
m = 0.06*(9.109e-31)
a = np.linspace(6e-9,12e-9,linspace_size)

InP_lc = (5.8697 - 6.0583)/(-0.4051)

InP_lc = round(InP_lc,2)
    
E_absorp = ((h*c)/(5e-6))/(eV)

m0 = (9.109e-31)
eff_me = ( (0.025*(1-InP_lc)) + (0.071*InP_lc) \
          - ((0.0163*InP_lc)*(1-InP_lc)) )*m0
    
E0e, roots1, roots2 = np.empty(len(a)), np.empty(len(a)), np.empty(len(a))

def Solve(x):
    return np.sqrt(((d*(np.pi**2))/4) - x**2) - x*(np.tan(x))
def Solve2(x):
    a = np.sqrt((((d*(np.pi**2))/4) - x**2))
    b = (x)*(np.cos(x)/np.sin(x))
    return a + b

Figure, (SubPlot1) = plt.subplots(1,1,constrained_layout=True)

c = 0
for i in a: 
    E0e_temp = (h**2)/(8*eff_me*(i**2))
    E0e[c] = E0e_temp
    d = (0.52*(eV))/E0e_temp
        
    zi = np.linspace(0,20,linspace_size)
    ye = np.sqrt(((d*(np.pi**2))/4) - zi**2)
    yh = (-zi)*(np.cos(zi)/np.sin(zi))
    y2 = zi*(np.tan(zi))

    ye_new = ye.astype(np.float)
    yh_new = yh.astype(np.float)
    y2new = y2.astype(np.float)

    ye_new[ye <= 0] = np.nan
    yh_new[yh <= 0] = np.nan
    y2new[y2 <= 0] = np.nan

    SubPlot1.plot(zi,ye_new)
    SubPlot1.plot(zi,y2new)
    SubPlot1.plot(zi,yh_new)
    
    if c >= 0 and c < (0.7*(+(len(a))/10)):
        EstInt1 = 1.2
        EstInt2 = np.array([2.1])
    elif c >= (0.7*(+(len(a))/10)) and c < (2.5*(+(len(a))/10)):
        EstInt1 = 1.3
        EstInt2 = np.array([2.3])
    elif c >= (2.5*(int(len(a))/10)) and c < (5*((len(a))/10)):
        EstInt1 = 1.3
        EstInt2 = np.array([2.4])
    elif c >= (5*((len(a))/10)) and c <= (7.5*((len(a))/10)):
        EstInt1 = 1.3
        EstInt2 = np.array([2.5])
    else:
        EstInt1 = 1.3
        EstInt2 = np.array([2.5])

    root1 = sc.fsolve(Solve,EstInt1)
    roots1[c] = root1
    root2 = sc.fsolve(Solve2,EstInt2)
    roots2[c] = root2
    c += 1
    
SubPlot1.set(xlabel = '$\u03BE$', \
             ylabel = 'q($\u03BE$), p($\u03BE$)', \
             title = 'q(\u03BE),p(\u03BE),r(\u03BE) vs \u03BE for QWIP', \
             ylim = ([0,5]), xlim = ([0,3]))
SubPlot1.grid(True)
plt.show()

E1e = ((4*(roots1**2))/(np.pi**2))*E0e
E2e = ((4*(roots2**2))/(np.pi**2))*E0e
'''
The 'delE' commented out here (line 112) can be used to calculate the 
absorption energies directly, although
it is more computationally resourceful just to calculate E2-E1 since these
have already been found.
'''
#delE = (2*((h_bar2))/(eff_me*(a**2)))*(roots2**2 - roots1**2)/eV
a = a/1e-9
E1e, E2e = E1e/eV, E2e/eV
delE = E2e-E1e

Figure, (ax) = plt.subplots(1,1,constrained_layout=True)
ax.plot(a,E1e)
ax.plot(a,E2e)
ax.plot(a,delE)
ax.axhline(y=E_absorp, linewidth='0.8', \
                 linestyle='--', color='orange')
ax.axvline(x=7.16, linewidth='0.8', \
                 linestyle='--', color='orange')
ax.grid(b=True, which='major', linestyle=':', linewidth='2')
ax.legend([r'$E_{1}$',r'$E_{2}$',\
           '$\u0394 E = E_{2} - E_{1}$',\
           '$E_{absorption}$ at a=7.16nm'], \
                loc = 'best', prop={'size': 20})
plt.minorticks_on()
ax.set_ylim(0,0.5)
ax.set_xlim(6,12)
ax.set(xlabel = 'a (nm)' ,\
       ylabel = '$E_{n}$ (eV)' ,)
ax.set(title = ''r'$E_{1}$  and $E_{2}$ vs Well Width wrt ' \
                r'the Conduction Band')
plt.show()

"""

23/02/2021

"""

linspace_size = 1000

a2 = 7.30e-9

E0 = ((h**2)/(8*eff_me*((a2**2))))
d = (0.52/E0)*eV

zi = np.linspace(0,20,linspace_size)
yelec = np.sqrt(((d*(np.pi**2))/4) - zi**2)
y1 = zi*(np.tan(zi))
y2 = ((-zi)*(np.cos(zi)/np.sin(zi)))

yelecnew = yelec.astype(np.float)
y1new = y1.astype(np.float)
y2new = y2.astype(np.float)

yelecnew[yelec <= 0] = np.nan
y1new[y1 <= 0] = np.nan
y2new[y2 <= 0] = np.nan

Figure, ax = plt.subplots(1,1,constrained_layout=True)
ax.plot(zi,yelecnew)
ax.plot(zi,y1new)
ax.plot(zi,y2new)
plt.ylim(0,5)
plt.xlim(0,5)
plt.show()

root1 = sc.fsolve(Solve,1.1)
root2 = sc.fsolve(Solve2,2)

x = np.linspace(-2*a2,2*a2,linspace_size)

psi1, psi2 = np.empty(len(x)), \
                    np.empty(len(x))

k1 = ((2*root1)/a2)
k1b = ((2*root2)/a2)
k2 = (2*(np.sqrt(((d*(np.pi**2))/4) - root1**2)))/a2
k2b = (2*(np.sqrt(((d*(np.pi**2))/4) - root2**2)))/a2

B = 1
A = 1

c = 0
for i in x:
    if x[c] < (-a2/2):
        psi1[c] = (B*(np.cos(root1))*(np.e**(k2*(a2/2))))*(np.e**(k2*i))
        psi2[c] = (-A*(np.sin(root2))*(np.e**(k2b*(a2/2))))*(np.e**(k2b*i))
    if x[c] > (-a2/2) and x[c] < (a2/2):
        psi1[c] = B*(np.cos( k1*i ))
        psi2[c] = A*(np.sin( k1b*i ))
    if x[c] > (a2/2):
        psi1[c] = (B*(np.cos(root1))*(np.e**(k2*(a2/2))))*(np.e**(-k2*i))
        psi2[c] = (A*(np.sin(root2))*(np.e**(k2b*(a2/2))))*(np.e**(-k2b*i))
    c += 1

Figure, ((ax,axx)) = plt.subplots(1,2,constrained_layout=True)
ax.plot(x/a2,psi1)
ax.plot(x/a2,psi2)
ax.grid(b=True, which='major', linestyle=':', linewidth='2')
ax.legend(['$\u03C8_{1}(x)$','$\u03C8_{2}(x)$'], \
                loc = 'best', prop={'size': 20})
ax.set_xlim(-2,2)
ax.set(xlabel = 'x/a',\
             ylabel = r'Normaslised Amplitude $(m^{(-\frac{1}{2})})$')
ax.set_title('The Eigenfunctions for the Bound Eigenstates ' \
                'wrt a=7.16nm', x=1)
axx.plot(x/a2,psi1**2)
axx.plot(x/a2,psi2**2)
axx.grid(b=True, which='major', linestyle=':', linewidth='2')
axx.set_xlim(-2,2)
axx.set(xlabel = 'x/a',\
             ylabel = r'Normaslised Probability Amplitude')
plt.show()

print("\n Take InGaAs Composition as In_{:.2f}Ga_{:.2f}As \n".format\
      ((1-InP_lc),InP_lc))

'''
#PREVIOUS GRAPH CODE

Figure, ((ax)) = plt.subplots(1,1,constrained_layout=True)
ax.plot(x/a2,psi1**2)
ax.plot(x/a2,psi2**2)
ax.grid(b=True, which='major', linestyle=':', linewidth='2')
ax.legend(['$\u03C8_{1}(x)$','$\u03C8_{2}(x)$'], \
                loc = 'best', prop={'size': 15})
ax.set_xlim(-2,2)
ax.set(xlabel = 'x/a',\
             ylabel = r'Normaslised Probability Amplitude',\
             title = 'The Eigenfunctions for the Bound Eigenstates ' \
                'wrt a=7.38nm')
plt.show()
'''

'''
#PLAYING WITH A


linspace_size = 1000

a2 = 1.3e-9

E0 = ((h**2)/(8*eff_me*((a2**2))))
d = (0.52/E0)*eV

zi = np.linspace(0,20,linspace_size)
yelectron = np.sqrt(((d*(np.pi**2))/4) - zi**2)
y1 = zi*(np.tan(zi))
y2 = ((-zi)*(np.cos(zi)/np.sin(zi)))

yelectronnew = yelectron.astype(np.float)
y1new = y1.astype(np.float)
y2new = y2.astype(np.float)

yelectronnew[yelectron <= 0] = np.nan
y1new[y1 <= 0] = np.nan
y2new[y2 <= 0] = np.nan

Figure, ax = plt.subplots(1,1,constrained_layout=True)
ax.plot(zi,yelectronnew)
ax.plot(zi,y1new)
ax.plot(zi,y2new)
plt.ylim(0,5)
plt.xlim(0,5)
plt.show()

root1 = sc.fsolve(Solve,0.35)

x = np.linspace(-3*a2,3*a2,linspace_size)

psi1, psi2 = np.empty(len(x)), \
                    np.empty(len(x))

k1 = ((2*root1)/a2)
k2 = (2*(np.sqrt(((d*(np.pi**2))/4) - root1**2)))/a2

#B = np.sqrt((2*k1*k2/(2*k1*(np.cos(k1*a2/2)**2)+k2*(np.sin(k1*a2)+k1*a2))))
B = 1
A = 1

c = 0
for i in x:
    if x[c] < (-a2/2):
        psi1[c] = (B*(np.cos(root1))*(np.e**(k2*(a2/2))))*(np.e**(k2*i))
    if x[c] > (-a2/2) and x[c] < (a2/2):
        psi1[c] = B*(np.cos( k1*i ))
    if x[c] > (a2/2):
        psi1[c] = (B*(np.cos(root1))*(np.e**(k2*(a2/2))))*(np.e**(-k2*i))
    c += 1

Figure, ((ax)) = plt.subplots(1,1,constrained_layout=False)
ax.plot(x/a2,psi1)
ax.grid(b=True, which='major', linestyle=':', linewidth='2')
ax.legend(['$\u03C8_{1}(x)$','$\u03C8_{2}(x)$'], \
                loc = 'best', prop={'size': 15})
#ax.set_ylim(0,10)
ax.set_xlim(-2,2)
ax.set(xlabel = 'x/a',\
             ylabel = r'Normaslised Amplitude $(m^{(-\frac{1}{2})})$',\
             title = 'The Eigenfunctions for the Bound Eigenstates ' \
                'wrt a=7.30nm')
plt.show()

Figure, ((ax)) = plt.subplots(1,1,constrained_layout=False)
ax.plot(x/a2,psi1**2)
ax.grid(b=True, which='major', linestyle=':', linewidth='2')
ax.legend(['$\u03C8_{1}(x)$','$\u03C8_{2}(x)$'], \
                loc = 'best', prop={'size': 15})
#ax.set_ylim(0,10)
ax.set_xlim(-3,3)
ax.set(xlabel = 'x/a',\
             ylabel = r'Normaslised Probability Amplitude',\
             title = 'The Eigenfunctions for the Bound Eigenstates ' \
                'wrt a=1.30nm')
plt.show()
'''