# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 12:15:11 2018

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
addnoise=0
RCseries=0
A= 5 #Amp

omega0=450*1e6*2*np.pi
#resonant frquency, where the difference in the two path is such that the 2 waves interfere destructively
tmax=1e-7
tstep=1e-10
t = np.arange(0,tmax ,tstep )
trange=int(tmax/tstep)
#Gamma0 calculation
Z0=50 +0j#Ohm
Cdev=32* 1e-18
Rc=1e3
#RC series

#RC parallel 
if RCseries==1:
    Z=Rc-1j/(omega0*Cdev)
else :   
 Z=Rc/(1+1j*omega0*Rc*Cdev) 
Gamma0=(Z-Z0)/(Z+Z0)
#plt.plot(linspace(400e6,500e6,1000),testsignfft)
#the phase offset is frequency dependent and it is annulled only at a particular frequency

lenght=10000


freq=[0]*lenght
Amplitude=[0]*lenght
Phase=[0]*lenght
offset=900e6 #at450MHz omega/offset=pi
#offset=offset/5#at450MHz omega/offset=5pi=pi
#I can already define phase and amplitude of cancellation signal that are not gonna change during the frequency scan
 #Bt is cancellation signal and is not affected by changement of capacitance of the DUT
phiav=[0]*lenght#int(tmax/100)
Z0=50 +0j#Ohm
Cdev=32* 1e-18
Rc=1e8
#simulate perfect cancellation
for i in range(0,lenght):
    
 omega= 400*1e6*2*np.pi +100*1e6*2*np.pi*i/lenght 
 Bt=(A*Gamma0)*exp(1j*(omega*t+ omega/offset))
 if addnoise==1:
  noise = np.random.normal(-0.001, 0.001, Bt.shape)
  Bt=Bt+noise
 #Bt=A*Gamma0/2*exp(1j*(omega*t+ omega/offset))
 #Bt=A*Gamma0*np.sin(omega*t+ omega/offset)
 #Define new reflection coefficient for different frequencies(Fixing Cdev)
# phi=1/np.arctan(-1/(omega*Cdev*Rc))
# ampZ=(Rc**2 +(omega*Cdev)**2)**0.5
# Z=ampZ* exp(1j*phi)
 #RC series
 #Z=Rc-1j/(omega*Cdev)
#RC parallel 
 if RCseries==1:
  Z=Rc-1j/(omega*Cdev)
 else :   
  Z=Rc/(1+1j*omega*Rc*Cdev) 

 Gamma=(Z-Z0)/(Z+Z0)

#Reflected signal
 #At=A*Gamma*np.sin(omega*t)
 At=A*Gamma*exp(1j*(omega*t))
 #At=A*Gamma/2 *exp(1j*(omega*t))
  #Bt= A*Gamma* (exp(1j*(omega*t+ np.pi)) - exp(-1j*(omega*t+ np.pi)))/(2j)

#cancellation condition
 freq[i]=400*1e6 +100*1e6*i/lenght
 cancelled= At+Bt
 signal=cancelled* (2**0.5)*exp(-1j*(omega-1)*t)
 I=np.real(signal)
 Q=np.imag(signal)
 
 Amp=(abs(I)**2 +abs(Q)**2)**0.5 
 phas=np.arctan2(Q,I) 

 Amplitude[i]=(Amp[trange])
 #Amplitude[i]=(Amp[trange])
 Phase[i]=phas[trange]
 Amplitude[i]=20*np.log10(Amplitude[i]/(A*abs(Gamma0))) 
 #Amplitude[i]=20*np.log10(Amplitude[i]/(A*abs(Gamma0))**2) 
 
 
Phase=np.unwrap(Phase)
plt.plot(freq,Amplitude)
plt.figure()
plt.plot(freq,Phase)
 
 #calculate Q factor
f=450e6
threshold=np.min(Amplitude)+3
k0=np.argmin(Amplitude)
k=k0
while Amplitude[k]<threshold:
 k=k-1
 print(k)
kmin=k
fminus= 400e6+ float(kmin)/lenght*1e8
k=k0
while Amplitude[k]<threshold:
 k=k+1
 kmax=k
fplus= 400e6+ float(kmax)/lenght*1e8
df=fplus-fminus
Q=f/df
print(Q)






#Ifft=np.fft.fft(signal)
#frequency=linspace(0,100e6,trange+1)
#plt.plot(frequency,Ifft)







#Now I can go back to the cancellation point  (or close to it) and simulate what happen for small capacitance variation
omega=omega0#+100000*(2*np.pi)

Clenght=10
ampC=[0]*Clenght
phiC=[0]*Clenght
deltaC=[0]*Clenght
Z2=[0]*Clenght
if RCseries==1:
    Z=Rc-1j/(omega0*Cdev)
else :   
 Z=Rc/(1+1j*omega0*Rc*Cdev) 
Gamma0=(Z-Z0)/(Z+Z0)

#Cdevalready defined
#Cdev=10*1e-18
for i in range(0,Clenght):
#deltaB due to non perfect cancellation--> can be related with Q
 #deltaB=float(A)/100*i
 deltaB=0
 Bt=(A*Gamma0+float(A)/100)*exp(1j*(omega*t+ omega/offset))
 if addnoise==1:
  noise = np.random.normal(-0.01, 0.01, Bt.shape)
  Bt=Bt+noise
#The signal change in amplitude and phase for small variation of the device capacitance, the cancellation signal Bt stays the same 
 deltaC[i]= -25e-18+float(i)/Clenght*50e-18
 C=Cdev+deltaC[i]
# phi2=1/np.arctan(-1/(omega*C*Rc))
# ampZ2=(Rc**2 +(omega*C)**2)**0.5
# Z2[i]=ampZ2* exp(1j*phi2)
 #Z2[i]=Rc-1j/(omega*Cdev)
#RC parallel
 if RCseries==1:
  Z2[i]=Rc-1j/(omega0*C)
 else :   
  Z2[i]=Rc/(1+1j*omega*Rc*C)  
 
  
 Gamma2=(Z2[i]-Z0)/(Z2[i]+Z0)
 At=Gamma2*A*(exp(1j*omega*t))
 cancelled= At+Bt

 signal=cancelled* (2**0.5)*exp(-1j*(omega-1)*t)
 I=np.real(signal)
 Q=np.imag(signal)
 Amp=(abs(I)**2 +abs(Q)**2)**0.5
 #dB convertion 
 phas=np.arctan2(Q,I) 
 ampC[i]=Amp[trange]
 #ampC[i]=Amp[trange]**2 

 phiC[i]=np.mean(phas)
 ampC[i]=20*np.log10(ampC[i]/(A*abs(Gamma0)))
# ampC[i]=20*np.log10(ampC[i]/(A*abs(Gamma0))**2)  
PhiC=np.unwrap(phiC) 

plt.plot(deltaC,ampC) 
plt.figure()
plt.plot(deltaC,PhiC)

#cancelledfft=np.fft.fft(cancelled)
#freq=linspace(0,100e6,trange+1)
#plt.plot(freq,cancelledfft)
#The broadening of the signal sent to the DUT in the f requency domain determines the width of the interfero line



#DeltaAmp vs Q , fixed delta C, on resonance
omega=omega0 #+100e3*2*np.pi
deltaC= 10e-18
C=Cdev+deltaC
Qlenght=20
Rc=10e7
if RCseries==1:
     Z=Rc-1j/(omega0*Cdev)
else :   
     Z=Rc/(1+1j*omega0*Rc*Cdev)  
Gamma0=(Z-Z0)/(Z+Z0)

Aref=[0]*lenght
Pref=[0]*lenght
freq=[0]*lenght

Ampmin=[0]*Qlenght
Phasemin=[0]*Qlenght
Qfactor=[0]*Qlenght
deltamp=[0]*Qlenght
deltaphi=[0]*Qlenght
A0=[0.]*Qlenght
P0=[0.]*Qlenght
addnoise=0
#First estimate the base signal than compare it to what happen when the dot capacitance change
for i in range(0,Qlenght): 
    
   for n in range (0,lenght):
    omega= 400*1e6*2*np.pi +100*1e6*2*np.pi*n/lenght
    #Calculate signal before of phase transition--> C=Cdev-->background depends on Q? NO
#    if i!=0:
    Bt=(A*Gamma0+float(A)/(100*float(i+1)))*exp(1j*(omega*t+ omega/offset))  
 
#    else:
#        Bt=A*Gamma0*exp(1j*(omega*t+ omega/offset))   
        #compare with no capacitance variations--> Calculate reference amplitude and phase
    if RCseries==1:
     Z2=Rc-1j/(omega*Cdev)
    else :   
     Z2=Rc/(1+1j*omega*Rc*Cdev)     
    
    Gamma2=(Z2-Z0)/(Z2+Z0)
    At=Gamma2*A*(exp(1j*omega*t))
    freq[n]=400*1e6 +100*1e6*n/lenght
    cancelled= At+Bt
    signal=cancelled* (2**0.5)*exp(-1j*(omega-1)*t)
    I=np.real(signal)
    Q=np.imag(signal)
 
    Amp=(abs(I)**2 +abs(Q)**2)**0.5 
    phas=np.arctan2(Q,I) 

    Aref[n]=(Amp[trange])
   # Aref[n]=(Amp[trange])**2
    Pref[n]=phas[trange]
    #Aref[n]=20*np.log10(Aref[n]/(A*abs(Gamma0))**2) 
    Aref[n]=20*np.log10(Aref[n]/(A*abs(Gamma0)))

    
    if addnoise==1:
     noise = np.random.normal(-0.0001, 0.0001, Bt.shape)
     Bt=Bt+noise
    
   #First cycle determines the accuracy of the frequency scan   
    freq[n]=400*1e6 +100*1e6*n/lenght
    #I take in account the Gamma modified by FIXED deltaC
    if RCseries==1:
     Z=Rc-1j/(omega*C)
    else :   
     Z=Rc/(1+1j*omega*Rc*C) 
#RC parallel 
 #Z=Rc/(1+1j*omega*Rc) 
    Gamma=(Z-Z0)/(Z+Z0)

#Reflected signal
 #At=A*Gamma*np.sin(omega*t)
    At=A*Gamma*exp(1j*(omega*t))
    cancelled= At+Bt
    signal=cancelled* (2**0.5)*exp(-1j*(omega-1)*t)
    I=np.real(signal)
    Q=np.imag(signal)
 
    Amp=(abs(I)**2 +abs(Q)**2)**0.5 
    phas=np.arctan2(Q,I) 
    Amplitude[n]=(Amp[trange])
    #Amplitude[n]=(Amp[trange])**2
    Phase[n]=phas[trange]
    Amplitude[n]=20*np.log10(Amplitude[n]/(A*abs(Gamma0))) 
    #Amplitude[n]=20*np.log10(Amplitude[n]/(A*abs(Gamma0))**2) 
    
   if addnoise==0 :
    f=450e6
    threshold=np.min(Amplitude)+3
    k0=np.argmin(Amplitude)
    k=k0
    while Amplitude[k]<threshold:
     k=k-1
    kmin=k
    print(kmin)
    fminus= 400e6+ float(kmin)/lenght*1e8
    k=k0
    while Amplitude[k]<threshold:
     k=k+1
    kmax=k
    print(kmax)
    fplus= 400e6+ float(kmax)/lenght*1e8
    df=fplus-fminus
    Qfactor[i]=float(f)/float(df)
   
   Ampmin[i]=Amplitude[k0]
   Phasemin[i]=Phase[k0]
   
   A0[i]=np.min(Aref)
   index=int(lenght/2)
   P0[i]=Pref[index]
 
   deltamp[i]=Ampmin[i]-A0[i]
   deltaphi[i]=Phasemin[i]-P0[i]
   #ampnoiz=
   #SNRamp[i]=deltaAmp[i]/ noiz  
plt.figure()   
plt.plot(Qfactor,deltamp,'o', label= '$\delta$A vs Q')
plt.figure()
plt.plot(Qfactor,deltaphi,'o', label= '$\delta (\phi)$ vs Q')  
plt.legend()
plt.title('')
plt.xlabel('$\delta$ P (dB)')
plt.ylabel('Q')
plt.savefig('dA vs Q')

#print(Q)