# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 16:38:02 2018

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
t = np.arange(0, 1e-7, 1e-10)

Z0=50 +0j#Ohm
Cdev=10* 1e-18
Rc=1000
phi=1/np.arctan(-1/(omega0*Cdev*Rc))
ampZ=(Rc**2 +(omega0*Cdev)**2)
Z=ampZ* exp(1j*phi)

#plt.plot(linspace(400e6,500e6,1000),testsignfft)
#the phase offset is frequency dependent and it is annulled only at a particular frequency
lenght=100
amp =[0]*lenght
freq=[0]*lenght
Amplitude=[0]*lenght
phase=[0]*lenght
offset=np.pi

for i in range(0,lenght):
    
 omega= 400*1e6*2*np.pi +100*1e6*2*np.pi*i/lenght 
 #At is cancellation signal and is not affected by changement of capacitance of the DUT
 At=A*(exp(1j*(omega*t+ offset*omega)) - exp(-1j*(omega*t+ offset*omega)))/(2j)
 Z0=50 +0j#Ohm
 Cdev=10* 1e-18
 Rc=1000
 phi=1/np.arctan(-1/(omega*Cdev*Rc))
 ampZ=(Rc**2 +(omega*Cdev)**2)
 Z=ampZ* exp(1j*phi)
 Gamma=(Z-Z0)/(Z+Z0)

#cancellation signal matched with reflected signal for all the frequencies, but the phase is not matched
 At=At*Gamma
 
  


#Bt= A*Gamma* (exp(1j*(omega*t+ np.pi)) - exp(-1j*(omega*t+ np.pi)))/(2j)
 Bt= A*abs(Gamma)* (exp(1j*omega*t) - exp(-1j*omega*t ))/(2j)
#Bt= A*Gamma* np.sin(omega*t +np.pi)
#cancellation condition
 freq[i]=400*1e6 +100*1e6*i/lenght
 cancelled= At-Bt
 
 cancelledfft=np.fft.fft(cancelled)
 plt.plot(cancelledfft) 
 

#it's ok Gamma is complex and alterate the phase as well Cancelled signal is already not 0???
 
 #Lock-in detection
 #in-phase component
 I= cancelled*  A*Gamma* (exp(1j*omega*t) - exp(-1j*omega*t ))/(2j)
 #quadrature component
 Q= cancelled*  A*Gamma* np.cos(omega*t)
 #IQ demodulation
 Amp=(abs(I)**2 +abs(Q)**2)**0.5
 phas= np.arctan(Q/I)
 #low pass filter to cut HF component
 Amplitude[i]= np.mean(Amp)
 phase[i]=phase[lenght-50:lenght]
#20 * np.log10(np.sqrt(2) * np.sqrt(x[-1]**2+y[-1]**2)/self._uhf.signal_output1_amplitude.get())  
 
 
 #amp[i]=max(signal)
 
 
 
 
 
 
plt.plot(freq,Amplitude)

