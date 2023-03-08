# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 12:10:09 2018

@author: AA255540
"""
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp, linspace, random
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

x = []
y = np.zeros(10)
Q=np.zeros(10)
#filename=S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\resume_lowQ.txt
foldername="Z:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\"
filename="snr.csv"
with open(filename, "w") as outfile:
    outfile.write('Q; amplitude peak;sigma peak; SNR')
    outfile.write("\n")

for i in range(0,10):
 if i==0:   
  Q[i]=18000    
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\phisigma_100points_6dBm_2018-07-25_13-46-18.txt'
  peakphi=0.0328

 if i==1:   
  Q[i]=2800    
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\phisigma_100points_7dBm_2018-07-25_12-36-32.txt' 
  peakphi=0.0129
 if i==2:   
  Q[i]=23000  
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\phisigma_100points_7dBm_2018-07-25_13-19-51.txt'
  peakphi=0.0324
 if i==3:   
  Q[i]=9200  
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\phisigma_100points_5dBm_2018-07-25_14-12-59.txt'
  peakphi=0.0169

 if i==4:   
  Q[i]=3100  
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\phisigma_100points_5dBm_2018-07-25_15-24-14.txt'
  peakphi=0.0063

 if i==5:   
  Q[i]=1750    
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\phisigma_100points_5dBm_2018-07-25_16-13-29.txt'
  peakphi=0.0046

 if i==6:   
  Q[i]=950    
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\phisigma_100points_2dBm_2018-07-25_03-02-16.txt'
  peakphi=0.0028
 if i==7:   
  Q[i]=1000    
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\phisigma_100points_3dBm_2018-07-25_01-45-16.txt'
  peakphi=0.0025
 if i==8:   
  Q[i]=1460    
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-24\\phisigma_100points_6dBm_2018-07-24_21-55-43.txt'
  peakphi=0.0058

 if i==9:   
  Q[i]=1700    
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-24\\phisigma_100points_5dBm_2018-07-24_23-12-16.txt'
  peakphi=0.0045


#name='resume_lowQ.txt'
 with open(folder+name,'r') as f:
  lines = f.readlines()
  x1 = np.array([line.split()[0] for line in lines],dtype='float64')
  y1 = abs(np.array([line.split()[1] for line in lines],dtype='float64'))
  #y2=  abs(np.array([line.split()[2] for line in lines],dtype='float64'))
  sigmapeak=y1[58:62]
  averagesigma=np.mean(sigmapeak)
  snr=peakphi/averagesigma
  y[i]=snr
  x=Q
 
 
 with open(filename, "a") as outfile:
    outfile.write('%.2f;%.4f;%.4f;%.4f;' % (float(Q[i]),peakphi, averagesigma, snr))
    outfile.write("\n")
fig = plt.figure()
plt.plot(Q, y,"o")
#plt.semilogx(t, np.sin(2*np.pi*t)
ax1 = fig.add_subplot(111)
 #ax1.set_xscale("log", nonposx='clip')
#ax1.set_title("Plot title")    
ax1.set_xlabel('Q')
#ax1.set_ylabel('$\sigma(\phi)$')
#ax1.set_ylabel('$\sigma(A)$')
#ax1.set_ylabel('$\phi_{max}$ / $\sigma_{peak}(\phi)$')
ax1.set_ylabel('$SNR(\phi)$')
ax1.plot(x,y, "o")#, label='the data')

leg = ax1.legend()

plt.savefig('SNR_phase_onthepeak.png')
plt.show()    