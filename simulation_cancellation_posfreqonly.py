# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 15:55:36 2018

@author: AA255540
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import exp, linspace, random
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import cmath 
import freqdemod
#amplitude in dBm
A= 5

omega0=450*1e6*2*np.pi
#resonant frquency, where the difference in the two path is such that the 2 waves interfere destructively
tmax=1e-7
tstep=1e-10
t = np.arange(0,tmax ,tstep )
trange=int(tmax/tstep)
Z0=50 +0j#Ohm
Cdev=10* 1e-18
Rc=1000
phi0=1/np.arctan(-1/(omega0*Cdev*Rc))
ampZ=(Rc**2 +(omega0*Cdev)**2)**0.5
Z=ampZ* exp(1j*phi0)
Gamma0=(Z-Z0)/(Z+Z0)
#plt.plot(linspace(400e6,500e6,1000),testsignfft)
#the phase offset is frequency dependent and it is annulled only at a particular frequency
lenght=1000



freq=[0]*lenght
Amplitude=[0]*lenght
Phase=[0]*lenght
offset=900e6 #at450MHz omega/offset=pi
#I can already define phase and amplitude of cancellation signal that are not gonna change during the frequency scan
 #Bt is cancellation signal and is not affected by changement of capacitance of the DUT
phiav=[0]*lenght#int(tmax/100)

for i in range(0,lenght):
    
 omega= 400*1e6*2*np.pi +100*1e6*2*np.pi*i/lenght 
 Bt=A*Gamma0*exp(1j*(omega*t+ omega/offset))
 #Bt=A*Gamma0*np.sin(omega*t+ omega/offset)
 #Define new reflection coefficient for different frequencies(Fixing Cdev)
 Z0=50 +0j#Ohm
 Cdev=10* 1e-18
 Rc=1000
 phi=1/np.arctan(-1/(omega*Cdev*Rc))
 ampZ=(Rc**2 +(omega*Cdev)**2)**0.5
 Z=ampZ* exp(1j*phi)
 Gamma=(Z-Z0)/(Z+Z0)

#Reflected signal
 #At=A*Gamma*np.sin(omega*t)
 At=A*Gamma*exp(1j*(omega*t))
  #Bt= A*Gamma* (exp(1j*(omega*t+ np.pi)) - exp(-1j*(omega*t+ np.pi)))/(2j)

#cancellation condition
 freq[i]=400*1e6 +100*1e6*i/lenght
 cancelled= At+Bt
# module=np.abs(cancelled)
# arg=np.angle(cancelled)
# 
# Amplitude[i]=module[trange] 
# Phase[i]=arg[100]
 
 #cancelledfft=np.fft.fft(cancelled)
 #plt.plot(cancelled) 
 #amp[i]=np.max(abs(cancelled))
 #plt.plot(freq,amp)
#Reference amplitude is AGamma0 
 #I= cancelled*  A*Gamma0* np.sin((omega-1)*t)
 #quadrature component
 #Q= cancelled*  A*Gamma0* np.cos((omega-1)*t)

# signal= cancelled* np.sin(omega*t)
# I=np.real(signal)
# Q=np.imag(signal)

# quadrature component
 I=cancelled* np.sin(omega*t)
 Q= cancelled* np.cos(omega*t)
 #IQ demodulation
 Amp=(abs(I)**2 +abs(Q)**2)**0.5
 #dB convertion 
 phas=np.arctan2(np.real(Q),np.real(I)) 
 #phas= np.arctan2(Q[trange],I[trange])
 #low pass filter to cut HF component

 Amplitude[i]=np.mean(Amp) 
 #20 * np.log10(np.sqrt(2) * np.sqrt(x[-1]**2+y[-1]**2)/self._uhf.signal_output1_amplitude.get())  
 Phase[i]=np.mean(phas)
 #Phase[i]=phas[400]
#it's ok Gamma is complex and alterate the phase as well Cancelled signal is already not 0???
  
 #Considering directly higher frequency component of the sine wave MAYBE I dont need demodulation
 #Lock-in detection
 #in-phase component   
plt.plot(freq,Amplitude)
plt.plot(freq,Phase)

Ifft=np.fft.fft(cancelled)
freq=linspace(0,100e6,trange+1)
plt.plot(freq,Ifft)



#Now I can go back to the cancellation point and simulate what happen for small capacitance variation
omega=omega0#+10000
#INTERESTING FLAT AT OMEGA0, linear at OMEGA0+10kHz
Clenght=10
ampC=[0]*Clenght
phiC=[0]*Clenght
deltaC=[0]*Clenght
Z2=[0]*Clenght
#Cdevalready defined
#Cdev=10*1e-18
for i in range(0,Clenght):
#deltaB due to non perfect cancellation--> can be related with Q
 #deltaB=float(A)/100*i
 deltaB=0
 Bt=A*Gamma0*exp(1j*(omega*t+ omega/offset))
#The signal change in amplitude and phase for small variation of the device capacitance, the cancellation signal Bt stays the same 
 deltaC[i]= i*1e-18
 C=Cdev+deltaC[i]
 phi2=1/np.arctan(-1/(omega*C*Rc))
 ampZ2=(Rc**2 +(omega*C)**2)**0.5
 #High frequency signal with + sign
 Z2[i]=ampZ2* exp(1j*phi2)
 Gamma2=(Z2[i]-Z0)/(Z2[i]+Z0)
 At=Gamma2*A*(exp(1j*omega*t))
 cancelled= At+Bt
 module=np.abs(cancelled)
 arg=np.angle(cancelled)
 ampC[i]=module[trange] 
 phiC[i]=arg[trange]
#lock in
# I= cancelled*  A*Gamma0* np.sin((omega-1)*t)
# #quadrature component
# Q= cancelled*  A*Gamma0* np.cos((omega-1)*t)
# #IQ demodulation
# Amplitude=(abs(I)**2 +abs(Q)**2)**0.5
# phas= np.arctan(Q/I)
 #low pass filter to cut HF component
 #Amplitude[i]= np.mean(Amp)
 #Here Amp and phase evaluated in a random point, no need to average
 #For the amplitude nothing change
 #For the phase the results depends on which time I check it, not just on how long
 #ampC[i]=Amplitude[trange] 
 #20 * np.log10(np.sqrt(2) * np.sqrt(x[-1]**2+y[-1]**2)/self._uhf.signal_output1_amplitude.get())  
 #Phase[i]=np.mean(phas)
 #phiC[i]=phas[trange]
 
plt.plot(deltaC,ampC) 
plt.plot(deltaC,phiC)

cancelledfft=np.fft.fft(cancelled)
freq=linspace(0,100e6,trange+1)
plt.plot(freq,cancelledfft)