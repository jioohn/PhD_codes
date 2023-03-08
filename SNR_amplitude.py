# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 12:47:55 2018

@author: AA255540
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import exp, linspace, random
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

x = []
y = []
Q=np.zeros(9)
sigm=np.zeros(9)
snr=np.zeros(9)
#filename=S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\resume_lowQ.txt
foldername="S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\"
filename="snr2.csv"
with open(foldername+filename, "w") as outfile:
    outfile.write('Q; amplitude peak;sigma peak; SNR')
    outfile.write("\n")

for i in range(0,9):
 if i==0:   
  Q[i]=18000    
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\phicorrect_100points_6dBm_2018-07-25_13-46-19.txt'
  peakphi=0.0328

 if i==1:   
  Q[i]=2800    
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\phicorrect_100points_7dBm_2018-07-25_12-36-33.txt' 
  peakphi=0.0129
 if i==2:   
  Q[i]=23000  
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\phicorrect_100points_7dBm_2018-07-25_13-19-52.txt'
  peakphi=0.0324
 if i==3:   
  Q[i]=9200  
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\phicorrect_100points_5dBm_2018-07-25_14-13-00.txt'
  peakphi=0.0169

 if i==4:   
  Q[i]=3100  
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\phicorrect_100points_5dBm_2018-07-25_15-24-15.txt'
  peakphi=0.0063

 if i==5:   
  Q[i]=1750    
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\phicorrect_100points_5dBm_2018-07-25_16-13-30.txt'
  peakphi=0.0046

 if i==6:   
  Q[i]=950    
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\correct_100points_2dBm_2018-07-25_03-02-17.txt'
  peakphi=0.0028
 if i==7:   
  Q[i]=1000    
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-25\\correct_100points_3dBm_2018-07-25_01-45-17.txt'
  peakphi=0.0025
 if i==8:   
  Q[i]=1460    
  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
  name='2018-07-24\\correct_100points_6dBm_2018-07-24_21-55-44.txt'
  peakphi=0.0058

# if i==9:   
#  Q[i]=1700    
#  folder='S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\data\\'
#  name='2018-07-24\\correct_100points_5dBm_2018-07-24_23-12-17.txt'
#  peakphi=0.0045


#name='resume_lowQ.txt'
 with open(folder+name,'r') as f:
  lines = f.readlines()
  x = np.array([line.split()[0] for line in lines],dtype='float64')
  y = np.array([line.split()[1] for line in lines],dtype='float64')
#  #y2=  abs(np.array([line.split()[2] for line in lines],dtype='float64'))
#  sigmapeak=y1[58:62]
#  averagesigma=np.mean(sigmapeak)
#  snr=peakphi/averagesigma
#  y[i]=snr
#  x=Q
 Q2=str(Q[i])
 def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))
 
 n = len(x)                          #the number of data

 mean = sum(x*y)/sum(y)                   #note this correction
 sigma =  np.sqrt(sum(y * (x - mean)**2) / sum(y))     #note this correction


 #GAUSSIAN FIT
 #popt,pcov = curve_fit(gaus,x,y,p0=[max(y),mean,sigma])
# residual=sum( (gaus(x,*popt) - y)**2)/n
# sigma= np.sqrt(residual)
# snr[i]=peakphi/sigma
# variance=residual - sum( (gaus(x,*popt) - y))**2    
# sigma= np.sqrt(variance)
# snr[i]=peakphi/sigma
# 
# plt.plot(x,y,'b+:',label='data')
# plt.plot(x,gaus(x,*popt),'ro:',label='fit')
# plt.legend()
# plt.title('Gaussian fit, Q='+Q2)
# plt.xlabel('Voltage (mV)')
# plt.ylabel('$\phi$ (mrad)')
# plt.savefig('gaussian_phasefit2_Q='+Q2+'.png')
# plt.show()
#LORENTZIAN FIT
 def lorentz(x,*p):
   return p[0] /((3.14159/2)*p[2])* 1/(1.0 + ((x  - p[1])/p[2])**2)
#def lorentz(x,*p):
#        return (1/(1+(x/p[2] - 1)**4*p[1]**2))*p[0]
    
 popt,pcov = curve_fit(lorentz,x,y,p0=[max(y),mean,sigma])
 residual=sum( (lorentz(x,*popt) - y)**2)/n
 #residual=sum( (gaus(x,*popt) - y)**2)/n
 sigma= np.sqrt(residual)
 snr[i]=peakphi/sigma

 plt.plot(x,y,'b+:',label='data')
 plt.plot(x,lorentz(x,*popt),'ro:',label='fit')
 plt.legend()
 plt.title('Lorentian fit, Q='+Q2)
 plt.xlabel('Voltage (mV)')
 plt.ylabel('\phi (rad)')
 plt.savefig('lorentian_phifit_Q='+Q2+'.png')
 plt.show()
 perr = np.sqrt(np.diag(pcov))
 width=2.355*popt[2]
 erwidth=2.355*perr[2]
 
 
 
 
 with open(foldername+filename, "a") as outfile:
    outfile.write('%.2f;%.6f;%.6f;%.4f;' % (float(Q[i]),peakphi, sigma, float(snr[i])))
    outfile.write("\n")
 sigm[i]=sigma

fig = plt.figure()
plt.plot(Q, snr,"o")
#plt.semilogx(t, np.sin(2*np.pi*t)
ax1 = fig.add_subplot(111)
 #ax1.set_xscale("log", nonposx='clip')
#ax1.set_title("Plot title")    
ax1.set_xlabel('Q')
#ax1.set_ylabel('$\sigma(\phi)$')
#ax1.set_ylabel('$\sigma(A)$')
#ax1.set_ylabel('$\phi_{max}$ / $\sigma_{peak}(\phi)$')
ax1.set_ylabel('$SNR (\phi)$')
#ax1.set_yscale("log", nonposy='clip')
ax1.plot(Q,snr, "o")#, label='the data')

leg = ax1.legend()

plt.savefig('Qvs_SNR_lorentzianfit_phi.png')
plt.show()