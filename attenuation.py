# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 15:19:10 2021

@author: AA255540
"""



kb=8.617333262145e-5#eV\K
h=6.582119569e-13#ev*s

# C=220e-12
# R=300e3

# R=1e6
# C=10e-9


# f0=1/(2*np.pi*R*C)
# print(f0)


f=5e9

T1=300
nb1=(h*f)/(np.exp(h*f/(kb*T1))-1)

T2=50
nb2=(h*f)/(np.exp(h*f/(kb*T2))-1)

T3=4
nb3=(h*f)/(np.exp(h*f/(kb*T3))-1)

T4=0.4
nb4=(h*f)/(np.exp(h*f/(kb*T4))-1)

# print(10*np.log10(nb1/nb2))

# print(10*np.log10(nb2/nb3))

# print(10*np.log10(nb3/nb4))


# print(10*np.log10(nb1/nb4))


A1=1.6
A2=100
A3=2



noise1=nb1/A1+(A1-1)*nb2/A1
noise2=noise1/A2+(A2-1)*nb3/A2

noise3=noise2/A3+(A3-1)*nb4/A3


noise3=noise3/kb

print(noise3)

kb=1.380649e-23 #J/K
e=1.602176634e-19#C
T=300
SS=kb*T/e*1e3*np.log(10)#mV




