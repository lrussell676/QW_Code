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

Task Five - 18/02/2021

"""

linspace_size = 10000
c = (3e8)
eV = (1.66e-19)
h = 6.62607015e-34
h_2 = (6.62607015e-34)**2
h_bar2 = ((6.62607015e-34)/(2*np.pi))**2
m = 0.06*(9.11e-31)
a = np.linspace(0.1e-9,10e-9,linspace_size)

InP_lc = (5.8697 - 6.0583)/(-0.4051)

print("\n Take InGaAs Composition as In_{:.4f}Ga_{:.4f}As \n".format\
      (InP_lc,(1-InP_lc)))

InP_lc = round(InP_lc)
    
E_absorp = ((h*c)/(5e-6))/(eV)

m0 = (9.11e-31)
eff_me = ( (0.025*(1-InP_lc)) + (0.071*InP_lc) \
          - ((0.0163*InP_lc)*(1-InP_lc)) )*m0
    
E0e, roots = np.empty(len(a)), np.empty(len(a))

def Solve(x):
    return np.sqrt(((d*(np.pi**2))/4) - x**2) - x*(np.tan(x))
def Solve2(x):
    cunt1 = np.sqrt((((d*(np.pi**2))/4) - x**2))
    cunt2 = (x)*(np.cos(x)/np.sin(x))
    return cunt1 + cunt2

c = 0
for i in a: 
    E0e_temp = (h**2)/(8*eff_me*(i**2))
    E0e[c] = E0e_temp
    d = (0.52*(eV))/E0e_temp
    if c > ((len(a))/10) and c <= (2*((len(a))/10)):
        EstInt = np.array([2.15])
    elif c> (2*((len(a))/10)) and c <= (3*((len(a))/10)):
        EstInt = np.array([2.15])
    elif c > (3*((len(a))/10)) and c <= (4*((len(a))/10)):
        EstInt = np.array([2.15])
    elif c > (4*((len(a))/10)) and c <= (5*((len(a))/10)):
        EstInt = np.array([2.15])
    elif c > (5*((len(a))/10)) and c <= (6*((len(a))/10)):
        EstInt = np.array([2.15])
    elif c > (6*((len(a))/10)) and c <= (7*((len(a))/10)):
        EstInt = np.array([2.15])
    elif c > (7*((len(a))/10)) and c <= (8*((len(a))/10)):
        EstInt = np.array([2.15])
    elif c > (8*((len(a))/10)) and c <= (9*((len(a))/10)):
        EstInt = np.array([2.16])
    else:
        EstInt = np.array([2.17])        
    
    root = sc.fsolve(Solve,EstInt)
    roots[c] = root
    c += 1

E1e = ((4*roots**2)/(np.pi**2))*E0e
a = a/1e-9
E0e, E1e = E0e/eV, E1e/eV
delE = E1e-E0e
Figure, ((ax,axx)) = plt.subplots(2,1,constrained_layout=False)
ax.plot(a,E0e)
ax.plot(a,E1e)
ax.plot(a,delE)
'''
ax.plot(a,delE)#ax.axhline(y=E_absorp, linewidth='0.8', \
                 linestyle='--', color='orange')
ax.axvline(x=6.16, linewidth='0.8', \
                 linestyle='--', color='orange')
ax.set(xlabel = 'a (nm)',\
             ylabel = '$E_{n}$ (eV)',\
             title = r'$E_{0}$  and $E_{1}$ vs Well Width wrt ' \
                r'the Conduction Band')
'''
ax.axhline(y=E_absorp, linewidth='0.8', \
                 linestyle='--', color='orange')
ax.axvline(x=7.30, linewidth='0.8', \
                 linestyle='--', color='orange')
ax.grid(b=True, which='major', linestyle=':', linewidth='2')
ax.legend([r'$E_{0}$',r'$E_{1}$',\
           '$\u0394 E = E_{1} - E_{0}$',\
           '$E_{absorption}$ at a=7.30nm'], \
                loc = 'best', prop={'size': 15})
axx.plot(a,E0e)
axx.plot(a,E1e)
axx.plot(a,delE)
axx.axhline(y=E_absorp, linewidth='0.8', \
                 linestyle='--', color='orange')
axx.axvline(x=7.30, linewidth='0.8', \
                 linestyle='--', color='orange')
axx.grid(b=True, which='major', linestyle=':', linewidth='2')
plt.minorticks_on()
axx.grid(b=True, which='minor', linestyle='-', linewidth='0.2')
ax.set_ylim(0,10)
ax.set_xlim(1,10)
axx.set_ylim(0,0.8)
axx.set_xlim(7,8)
Figure.text(0.5, 0.94, ''r'$E_{0}$  and $E_{1}$ vs Well Width wrt ' \
                r'the Conduction Band', ha='center', fontsize=16)
Figure.text(0.5, 0.035, 'a (nm)', ha='center', fontsize=15)
Figure.text(0.055, 0.5, '$E_{n}$ (eV)', va='center', rotation='vertical', \
            fontsize=15)
plt.show()

"""

Task Five C - 23/02/2021

"""

linspace_size = 1000

a = 7.30e-9

E0 = ((h**2)/(8*eff_me*((a**2))))
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

root1 = sc.fsolve(Solve,1.1)
root2 = sc.fsolve(Solve2,2)

E1 = ((root1/4)**2)
E2 = ((root2/4)**2)

x = np.linspace(-2*a,2*a,linspace_size)

psi1, psi2 = np.empty(len(x)), \
                    np.empty(len(x))

k1 = ((2*root1)/a)
k1b = ((2*root2)/a)
k2 = (2*(np.sqrt(((d*(np.pi**2))/4) - root1**2)))/a
k2b = (2*(np.sqrt(((d*(np.pi**2))/4) - root2**2)))/a

#B = np.sqrt((2*k1*k2/(2*k1*(np.cos(k1*a/2)**2)+k2*(np.sin(k1*a)+k1*a))))
B = 1
A = 1

c = 0
for i in x:
    if x[c] < (-a/2):
        psi1[c] = (B*(np.cos(root1))*(np.e**(k2*(a/2))))*(np.e**(k2*i))
        psi2[c] = (-A*(np.sin(root2))*(np.e**(k2b*(a/2))))*(np.e**(k2b*i))
    if x[c] > (-a/2) and x[c] < (a/2):
        psi1[c] = B*(np.cos( k1*i ))
        psi2[c] = A*(np.sin( k1b*i ))
    if x[c] > (a/2):
        psi1[c] = (B*(np.cos(root1))*(np.e**(k2*(a/2))))*(np.e**(-k2*i))
        psi2[c] = (A*(np.sin(root2))*(np.e**(k2b*(a/2))))*(np.e**(-k2b*i))
    c += 1

Figure, ((ax)) = plt.subplots(1,1,constrained_layout=False)
ax.plot(x/a,psi1)
ax.plot(x/a,psi2)
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
ax.plot(x/a,psi1**2)
ax.plot(x/a,psi2**2)
ax.grid(b=True, which='major', linestyle=':', linewidth='2')
ax.legend(['$\u03C8_{1}(x)$','$\u03C8_{2}(x)$'], \
                loc = 'best', prop={'size': 15})
#ax.set_ylim(0,10)
ax.set_xlim(-2,2)
ax.set(xlabel = 'x/a',\
             ylabel = r'Normaslised Probability Amplitude',\
             title = 'The Eigenfunctions for the Bound Eigenstates ' \
                'wrt a=7.30nm')
plt.show()


#PLAYING WITH A


linspace_size = 1000

a = 1.30e-9

E0 = ((h**2)/(8*eff_me*((a**2))))
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

root1 = sc.fsolve(Solve,0.3)

E1 = ((root1/4)**2)

x = np.linspace(-3*a,3*a,linspace_size)

psi1, psi2 = np.empty(len(x)), \
                    np.empty(len(x))

k1 = ((2*root1)/a)
k2 = (2*(np.sqrt(((d*(np.pi**2))/4) - root1**2)))/a

#B = np.sqrt((2*k1*k2/(2*k1*(np.cos(k1*a/2)**2)+k2*(np.sin(k1*a)+k1*a))))
B = 1
A = 1

c = 0
for i in x:
    if x[c] < (-a/2):
        psi1[c] = (B*(np.cos(root1))*(np.e**(k2*(a/2))))*(np.e**(k2*i))
    if x[c] > (-a/2) and x[c] < (a/2):
        psi1[c] = B*(np.cos( k1*i ))
    if x[c] > (a/2):
        psi1[c] = (B*(np.cos(root1))*(np.e**(k2*(a/2))))*(np.e**(-k2*i))
    c += 1

Figure, ((ax)) = plt.subplots(1,1,constrained_layout=False)
ax.plot(x/a,psi1)
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
ax.plot(x/a,psi1**2)
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
