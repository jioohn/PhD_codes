# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 17:21:57 2021

@author: AA255540
"""

import csv
import matplotlib.pyplot as plot

import numpy as np

from scipy.optimize import curve_fit
from scipy import optimize
from scipy import asarray as ar,exp
from scipy.stats import norm
from qcodes import validators as vals
import scipy
from scipy.optimize import curve_fit
from scipy import fftpack


kb=8.617333262145e-5
bohr=5.7883818012e-5
h = 4.135667662e-15 # in eV.s




####
# 1/T= omega_L**2* coth(hbar * omega_L)/(2*kb*T))
####

g=1.5
B=np.linspace(0.5,1,801)
fLarmor=g*bohr*B/h

factor=1e22


T=0.1
T1_2=fLarmor**-2*np.tanh(h*fLarmor/(2*kb*T))*factor


T=0.1
T1=fLarmor**-3*np.tanh(h*fLarmor/(2*kb*T))*factor 


plt.figure()
plt.semilogy(B,T1_2)


##linear fit

def fit_B(x,*p):
    
    
    return p[1]+p[2]*B**p[0]
    


popt=[]
pcov=[]

popt,pcov = curve_fit(fit_B,B,T1,p0=[-2,0,0],maxfev=10000)

plt.figure()
plt.plot(B,T1,label='simulation')
plt.plot(B,fit_B(B,*popt),label=r'slope ='+str(np.round(popt[0],2)))
plt.legend()
plt.savefig('FIT_T1vsBfield_simulation')




#### Temperature dependence
g=1.5
B=np.linspace(0.5,1,101)
fLarmor=g*bohr*B/h

factor=1e22


T=0.1
T1_2=fLarmor**-2*np.tanh(h*fLarmor/(2*kb*T))*factor


alpha=0.5e-10
beta=1
T=0.1
# T1=fLarmor**-2*np.tanh(h*fLarmor/(2*kb*T))*factor

Gamma1=alpha*fLarmor**3*1/np.tanh(h*fLarmor/(2*kb*T))+ beta*fLarmor**2*1/np.tanh(h*fLarmor/(2*kb*T))


T1=1/Gamma1

plt.figure()
plt.semilogy(B,T1_2)


##polynomial

def fit_B(x,*p):
    
    
    return p[1]+p[2]*B**p[0]
    


popt=[]
pcov=[]

popt,pcov = curve_fit(fit_B,B,T1,p0=[-2.5,1.50e-6,10e-6],maxfev=100000)

plt.figure()
plt.plot(B,T1,'+',label='simulation')
plt.plot(B,fit_B(B,*popt),label=r'slope ='+str(np.round(popt[0],2)))
plt.legend()
plt.savefig('FIT_BfieldvsT1_simulation')





def good_fit_B(B,*p):
    
    fLarmor=g*bohr*B
    Gamma2D= p[1]+p[0]*fLarmor**2*1/np.tanh(h*fLarmor/(2*kb*T))
    Gamma3D=p[2]+p[3]*fLarmor**3*1/np.tanh(h*fLarmor/(2*kb*T))
    
    return 1/(Gamma2D+Gamma3D)
    


#####bosonic DOS

T=0.1
g=1.5
B=1
fLarmor=g*bohr*B/h



# n1=1/(np.exp(-h*fLarmor/kb*T)-1)

# T=0.5


# n2=1/(np.exp(-h*fLarmor/kb*T)-1)

n1=1/(np.exp(h*fLarmor/kb*T)-1)

T=0.5


n2=1/(np.exp(h*fLarmor/kb*T)-1)

print(n2/n1)

####T dependence

g=1.5
B=1
fLarmor=g*bohr*B/h

T=np.linspace(0,1,1000)
n1=1/(np.exp(-h*fLarmor/kb*T)-1)
# n1=1/(np.exp(-h*fLarmor/kb*T)-1)
T1=1/n1




plt.figure()
plt.plot(T,T1)
# plt.semilogy(T,T1)
#####################################




#Loss paper 2phonons at 
hbar=h/(2*np.pi)
c=4e3#m/s in In-As
c=3e3#m/s in GA-As
# In silicon bulk along 00
c_t=3.47e3
c_l=4.49e3


#found in https://periodictable.com/Properties/A/SoundSpeed.al.html
c=3e3

lamda=37e-9

Eph=hbar*c/(lamda)

Tcrit=hbar*c/(lamda*kb*2)


print(Tcrit)