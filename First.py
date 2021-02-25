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
#import PH380Functions as PH380
import scipy.optimize as sc
#from scipy.optimize import fsolve # Finding Zeros of Functions section


"""

Task Two - 21/01/2021

"""
linspace_size = 10000

xArray = np.linspace(0,8,linspace_size)

test_zi = xArray

q_zi = np.sqrt( 16 - test_zi**2 )
p_zi = test_zi*(np.tan(test_zi))
p_zi[:-1][np.diff(p_zi)<0]=np.nan # the negative 1 here corrects array missize
r_zi = (-test_zi)*(np.cos(test_zi)/np.sin(test_zi))
r_zi[:-1][np.diff(r_zi)<0]=np.nan # the negative 1 here corrects array missize

Figure, Plot1 = plt.subplots(1,1)
Figure.patch.set_facecolor('xkcd:mint green')
Figure.suptitle("NumSolver for Eigenvalues")
Plot1.set(ylabel = "E", xlabel = "Eigenvalues")
Plot1.plot(xArray,q_zi)
Plot1.plot(xArray,p_zi)
Plot1.plot(xArray,r_zi)
Plot1.grid(b=True, which='major', linestyle=':', linewidth='2')

plt.ylim(0, 5)
plt.show()

EstimatedIntersectionP = np.array([1.3,3.8])
EstimatedIntersectionR = np.array([2.5,3.2])

def PToSolve(x): 
    return  (np.sqrt( 16 - x**2 )) - (x*(np.tan(x)))
def RToSolve(x): 
    return  (np.sqrt( 16 - x**2 )) - ((-x)*(np.cos(x)/np.sin(x)))  
SolutionsP = sc.fsolve(PToSolve,EstimatedIntersectionP)
SolutionsR = sc.fsolve(RToSolve,EstimatedIntersectionR)
print(f"\n \n The solutions are: {SolutionsP[0]} and {SolutionsP[1]}")
print(f"The solutions are: {SolutionsR[0]} and {SolutionsR[1]} \n \n")


"""

Task Three - 25/01/2021, 28/01/2021

"""

h_2 = (6.62607015e-34)**2
h_bar2 = ((6.62607015e-34)/(2*np.pi))**2
m = 0.06*(9.11e-31)

a = np.linspace(0.1e-9,10e-9,linspace_size)

E_0 = ((h_2/(8*m*(a**2))))/(1.6e-19)


Figure, Plot3 = plt.subplots(1,1)
Figure.patch.set_facecolor('xkcd:mint green')
Figure.suptitle("$E_{0}$ vs Well Width")
Plot3.set(ylabel = "$E_{0}$/eV", xlabel = "a (m)")
Plot3.plot(a,E_0)
Plot3.grid(b=True, which='major', linestyle=':', linewidth='2')

plt.ylim(0, 50)
plt.show()

d = np.array([0.1,1,10,100])

zi = np.linspace(0,19,linspace_size)
q_zi = np.empty((4,linspace_size))

for i in range(len(d)):
    q_zi[i,:] = np.sqrt((np.pi**2)*(d[i]/4) - zi**2)

p_zi = test_zi*(np.tan(zi))
p_zi[:-1][np.diff(p_zi)<0]=np.nan # the negative 1 here corrects array missize
r_zi = (-zi)*(np.cos(zi)/np.sin(zi))
r_zi[:-1][np.diff(r_zi)<0]=np.nan # the negative 1 here corrects array missize

Figure, ((Plot4,p4b),(p4c,p4d)) = plt.subplots(2,2,constrained_layout=True)
Figure.patch.set_facecolor('xkcd:mint green')
Figure.suptitle("E_n vs Well Depth (d)")
Plot4.set(ylabel = "E_0/eV", xlabel = "EigenValues at d (relative) = 0.1")
p4b.set(ylabel = "E_0/eV", xlabel = "EigenValues at d (relative) = 1")
p4c.set(ylabel = "E_0/eV", xlabel = "EigenValues at d (relative) = 10")
p4d.set(ylabel = "E_0/eV", xlabel = "EigenValues at d (relative) = 100")
Plot4.plot(zi,p_zi)
p4b.plot(zi,p_zi)
p4c.plot(zi,p_zi)
p4d.plot(zi,p_zi)
Plot4.plot(zi,r_zi)
p4b.plot(zi,r_zi)
p4c.plot(zi,r_zi)
p4d.plot(zi,r_zi)
Plot4.plot(zi,q_zi[0,:])
p4b.plot(zi,q_zi[1,:])
p4c.plot(zi,q_zi[2,:])
p4d.plot(zi,q_zi[3,:])
Plot4.grid(b=True, which='major', linestyle=':', linewidth='2')
p4b.grid(b=True, which='major', linestyle=':', linewidth='2')
p4c.grid(b=True, which='major', linestyle=':', linewidth='2')
p4d.grid(b=True, which='major', linestyle=':', linewidth='2')
Plot4.set_ylim(0, 0.5)
Plot4.set_xlim(0, 0.5)
p4b.set_ylim(0, 2)
p4b.set_xlim(0, 2)
p4c.set_ylim(0, 6)
p4c.set_xlim(0, 6)
p4d.set_ylim(0, 17.5)
plt.show()

def PToSolve3a(x): 
    return  (np.sqrt( (np.pi**2)*(d[0]/4) - x**2 )) - \
        (x*(np.tan(x)))
        
def RToSolve3a(x): 
    return  (np.sqrt( (np.pi**2)*(d[0]/4) - x**2 )) - \
        ((-x)*(np.cos(x)/np.sin(x)))
        
def PToSolve3b(x): 
    return  (np.sqrt( (np.pi**2)*(d[1]/4) - x**2 )) - \
        (x*(np.tan(x)))
        
def RToSolve3b(x): 
    return  (np.sqrt( (np.pi**2)*(d[1]/4) - x**2 )) - \
        ((-x)*(np.cos(x)/np.sin(x)))
        
def PToSolve3c(x): 
    return  (np.sqrt( (np.pi**2)*(d[2]/4) - x**2 )) - \
        (x*(np.tan(x)))
        
def RToSolve3c(x): 
    return  (np.sqrt( (np.pi**2)*(d[2]/4) - x**2 )) - \
        ((-x)*(np.cos(x)/np.sin(x)))
        
def PToSolve3d(x): 
    return  (np.sqrt( (np.pi**2)*(d[3]/4) - x**2 )) - \
        (x*(np.tan(x)))
        
def RToSolve3d(x): 
    return  (np.sqrt( (np.pi**2)*(d[3]/4) - x**2 )) - \
        ((-x)*(np.cos(x)/np.sin(x)))
        
EstInPa = (0.44)
EstInPb1 = (1.1)
EstInPc1 = (1.4)
EstInPc2 = (2.7)
EstInPc3 = (4.2)
EstInPc4 = (4.8)
EstInPd1 = (1)
EstInPd2 = (2)
EstInPd3 = (4)
EstInPd4 = (6)
EstInPd5 = (7)
EstInPd6 = (8)
EstInPd7 = (10.5)
EstInPd8 = (11)
EstInPd9 = (13)
EstInPd10 = (14)
EstInPd11 = (15)

SolutionPa = sc.fsolve(PToSolve3a,EstInPa)
SolutionPb1 = sc.fsolve(PToSolve3b,EstInPb1)
SolutionPc1 = sc.fsolve(PToSolve3c,EstInPc1)
SolutionPc2 = sc.fsolve(RToSolve3d,EstInPc2)
SolutionPc3 = sc.fsolve(PToSolve3c,EstInPc3)
SolutionPc4 = sc.fsolve(RToSolve3d,EstInPc4)
SolutionPd1 = sc.fsolve(PToSolve3d,EstInPd1)
SolutionPd2 = sc.fsolve(RToSolve3d,EstInPd2)
SolutionPd3 = sc.fsolve(PToSolve3d,EstInPd3)
SolutionPd4 = sc.fsolve(RToSolve3d,EstInPd4)
SolutionPd5 = sc.fsolve(PToSolve3d,EstInPd5)
SolutionPd6 = sc.fsolve(RToSolve3d,EstInPd6)
SolutionPd7 = sc.fsolve(PToSolve3d,EstInPd7)
SolutionPd8 = sc.fsolve(RToSolve3d,EstInPd8)
SolutionPd9 = sc.fsolve(PToSolve3d,EstInPd9)
SolutionPd10 = sc.fsolve(RToSolve3d,EstInPd10)
SolutionPd11 = sc.fsolve(PToSolve3d,EstInPd11)

print(f"\n \n The solutions A1 are: {SolutionPa[0]}")
print(f"\n The solutions B1 are: {SolutionPb1[0]}")
print(f"\n The solutions C1 are: {SolutionPc1[0]}")
print(f"\n The solutions C2 are: {SolutionPc2[0]}")
print(f"\n The solutions C3 are: {SolutionPc3[0]}")
print(f"\n The solutions C4 are: {SolutionPc4[0]}")
print(f"\n The solutions D1 are: {SolutionPd1[0]}")
print(f"\n The solutions D2 are: {SolutionPd2[0]}")
print(f"\n The solutions D3 are: {SolutionPd3[0]}")
print(f"\n The solutions D4 are: {SolutionPd4[0]}")
print(f"\n The solutions D5 are: {SolutionPd5[0]}")
print(f"\n The solutions D6 are: {SolutionPd6[0]}")
print(f"\n The solutions D7 are: {SolutionPd7[0]}")
print(f"\n The solutions D8 are: {SolutionPd8[0]}")
print(f"\n The solutions D9 are: {SolutionPd9[0]}")
print(f"\n The solutions D10 are: {SolutionPd10[0]}")
print(f"\n The solutions D11 are: {SolutionPd11[0]} \n \n")

sols = np.array([SolutionPa[0],SolutionPb1[0],SolutionPc1[0],SolutionPc2[0], \
                 SolutionPc3[0],SolutionPc4[0],SolutionPd1[0],SolutionPd2[0], \
                 SolutionPd3[0],SolutionPd4[0],SolutionPd5[0],SolutionPd6[0], \
                 SolutionPd7[0],SolutionPd8[0],SolutionPd9[0], \
                     SolutionPd10[0],SolutionPd11[0]])
    
E_E0 = (4*(sols**2)/(np.pi**2))
d_fake = np.array([0.1,1,10,10,10,10,100,100,100,100, \
                   100,100,100,100,100,100,100])
d_log_fake = np.log10(d_fake)

Figure, Plot3b = plt.subplots(1,1)
Figure.patch.set_facecolor('xkcd:mint green')
Figure.suptitle("Ratios")
Plot3b.set(ylabel = "$E_{n}$/$E_{0}$", xlabel = "Relative (log(10)) Depth")
Plot3b.plot(d_log_fake,E_E0, ls='', marker='o')
Plot3b.grid(b=True, which='major', linestyle=':', linewidth='2')
plt.show()


"""

Task Four - 04/02/2021, 11/02/2021

"""

al_ratio = np.linspace(0, 0.4, 100000)
E_g = (1.424 + 1.247*al_ratio)

Figure, Plot4 = plt.subplots(1,1)
Figure.patch.set_facecolor('xkcd:mint green')
Figure.suptitle("Al Ratio Plot")
Plot4.set(ylabel = "$E_{g}$ (eV)", xlabel = "AlGa Fraction (x)")
Plot4.plot(al_ratio,E_g)
Plot4.plot(al_ratio[83514],E_g[83514], marker='x', markersize=15)
#Plot4.plot(al_ratio[20128],E_g[20128], marker='x', markersize=15)
Plot4.grid(b=True, which='major', linestyle=':', linewidth='2')
Plot4.legend(['$E_{g}$','X Marker Point of E_pump','X Marker Point of E_emission'],loc='best')
plt.show()

h = (6.62607015e-34)
E_pump = (h* ((3e8)/(675e-9)) )/(1.6e-19)
E_emission = (h* ((3e8)/(815e-9)) )/(1.6e-19)

print(f"\n \n E_pump is {E_g[83514]}eV at AlGa Fraction {al_ratio[83514]}")

#print(f"\n E_emission is {E_g[20128]}eV at AlGa Fraction {al_ratio[20128]}")

print("\n Take 'A' Barrier Composition as Al_{0.330}Ga_{0.670}As")

print("\n Take 'B' Well Composition as GaAs")

Figure, Plot5 = plt.subplots(1,1)
Figure.patch.set_facecolor('xkcd:mint green')
Figure.suptitle("$E_{0}$ vs Well Width")
Plot5.set(ylabel = "$E_{0}$/eV", xlabel = "a (m)")
Plot5.plot(a,E_0)
Plot5.grid(b=True, which='major', linestyle=':', linewidth='2')


me = (9.11e-31)
mhh = 0.45*me
a = np.linspace(0.1e-9,7e-9,linspace_size)

E_0e = ((h_2/(8*me*(a**2))))/(1.6e-19)
E_0hh = ((h_2/(8*mhh*(a**2))))/(1.6e-19)

Figure, Plot6 = plt.subplots(1,1)
Figure.patch.set_facecolor('xkcd:mint green')
Figure.suptitle("$E_{0}$ vs Well Width")
Plot6.set(ylabel = "$E_{0}$/eV", xlabel = "a (nm)")
Plot6.plot(a*(1e9),E_0e, label = '$E_{e(0)}$ =' r'$\frac{h^2}{8m_{0}a^2}$')
Plot6.plot(a*(1e9),E_0hh, label = '$E_{hh(0)}$ =' r'$\frac{h^2}{8m_{hh}a^2}$')
Plot6.legend(loc = 'best')
Plot6.grid(b=True, which='major', linestyle=':', linewidth='2')
plt.minorticks_on()
Plot6.grid(b=True, which='minor', linestyle='-', linewidth='0.2')

plt.ylim(0, 0.5)
plt.show()

anm = np.sqrt( ((h_2*(1-0.46))/(8*me*((E_pump-E_emission)*(1.6e-19)) )))

print(anm)