# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 13:53:50 2021

@author: AA255540
"""


#T2=2*np.sqrt(np.log(2))/(np.pi*20e6)
epsilon0=8.854187e-12
epsilonr=3.9
A=75*40*1e-18
d=6e-9
Cox=epsilon0*epsilonr*A/d

hbar=6.626e-34
me=9.109e-31


#estimate area from DeltaVg=Ec/alpha
DeltaVG=0.0085#V
A=1.6e-19*6e-9/(epsilon0*epsilonr)/DeltaVG

L=np.sqrt(A)





A=np.pi*hbar**2/(me*300e-6*1.6e-19)
r=np.sqrt(A/np.pi)

#for 2Ddot 
L=60e-9
deltaN=np.pi*hbar**2/(me*L**2)/1.6e-19*1e3
print('spacing='+str(deltaN)+'meV')

Ten=kb*4.2*1e3
print('Ten'+str(Ten)+'meV')


#for Weyl forula
#] H.P. Baltes, E.R. Hilf, Spectra of Finite Systems: A Reviewof Weyl’s Problem, Bibliographisches Institut, Zurich, 1976
# R=20e-9
L=30e-9
# A=np.pi*R**2
hbar=1.054e-34
me=9.109e-31
me=0.32*me
g=4
deltaN=2*np.pi*hbar**2/(g*me*L**2)*1e3/1.6e-19
# deltaN=2*np.pi*hbar**2/(g*me*A)*1e3/1.6e-19
print('spacing='+str(deltaN)+'meV')




#] Kouwenhoven 3D dot
L=50e-9
hbar=1.054e-34
me=9.109e-31
g=4#spin degenereacy
N=1
deltaN=(1/(3*np.pi**2*N))**0.333*(np.pi**2*hbar**2)/(g*me*L**2)*1e3/1.6e-19
print('spacing_3Ddot='+str(deltaN)+'meV')





