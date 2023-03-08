#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 15:06:34 2023

@author: StefaniaParigi
"""


B=0.9

T=np.linspace(0.08,0.45,101)
T1=110*B**-3
fL=20e9

nb=1/(np.exp(h*fL/(kb*T)) -1)

plt.rcParams['font.size'] = '22'

B=0.9
g=1.5
fL=g*bohr*B/h
nb=1/(np.exp(h*fL/(kb*T)) -1)

plt.figure()
# plt.plot(T,nb)
plt.semilogy(T,nb**-1,label='B=1 T')

B=1.5
g=1.5
fL=g*bohr*B/h
nb=1/(np.exp(h*fL/(kb*T)) -1)

plt.semilogy(T,nb**-1,label='B=1.5 T')


plt.xlabel('T (K)')#,fontsize=18)
plt.ylabel('$n^{-1}_B$')
 
plt.savefig('boson_density_inv_log.png')


