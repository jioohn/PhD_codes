# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 16:51:12 2021

@author: AA255540
"""

e=1.6e-19#C
d=6e-9#m
epsr=3.9
alpha=0.1#eV/V
eps0=8.8541878128e-12#F/M
Ec=100e-3#{eV}
A=e*alpha/Ec*d/(eps0*epsr)

print(A**0.5)

e=1.6e-19#C
d=6e-9#m
epsr=3.9
alpha=0.27#eV/V
eps0=8.8541878128e-12#F/M
# Ec=20e-3*0.27#{eV}
Ec=6e-3
A=e*alpha/Ec*d/(eps0*epsr)

#diameter assuming planar circular dot
D=(A/np.pi)**0.5*2
print(D)

D=2*(A**0.5)/np.pi


#dot level spacing
h=4.135667696e-15#ev/s
hbar=h/np.pi
me=9.1e-31#Kg
h=6.62607015e-34#J/s
hbar=h/np.pi



# D=4*0.01*d/(eps0*epsr*np.pi)

delta=2*np.pi*hbar**2/(4*0.32*me* 50*70e-18)
