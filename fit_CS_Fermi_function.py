# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 16:13:15 2021

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

# 
# https://apps.automeris.io/wpd/
#calibrationFilePath = r'Z:\132.05.01-QuantumSilicon\Triton experiments\data\LineClibration\2019-12-06_09-53-33_LineCalibrationAt600mTtweaked.csv'
folder=r'C:\Users\AA255540\Desktop\Tritonito2021_alldata\Feb2021_6G22_2_die103\exploration\G3=-772.429_g2=-997.286\Default Dataset (13).csv'


calibrationFilePath=folder
calibs = []
with open(calibrationFilePath) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=';')
        line_count = 0
        for row in csv_reader:
                calibs.append(row)
                
                

x=[]
y=[]
for i in range (0,len(calibs)):
    x.append(float(calibs[i][0].replace(',','.')))
    y.append(float(calibs[i][1].replace(',','.')))
    
# plot.figure()
# plot.plot(x, y,'+',label='G1') 
# # plot.title('Is vs $V_{G}$')
# plot.xlabel('$V_{G3}$ (V)')
# plot.ylabel('$\phi$ (deg)')




# phase=scipy.signal.savgol_filter(y,11,5)
####fit

ph_ma = np.max(phase)
ph_mi = np.min(phase)
Delta = ph_mi - ph_ma 

y = phase/Delta-np.min(phase)/Delta
y=-y

plot.figure()
plot.plot(x, y,'+',label='data') 
# plot.title('Is vs $V_{G}$')
plot.xlabel('$V_{G3}$ (V)')
plot.ylabel('$\phi$ (deg)')
plot.legend()

# y=-y

T=0.44


# def tauvsfield_4(x,*p):
#  Ez=1.98*bohr*x
#  # gamma_down=p[0]*exp(-p[2]/(kb*T))*1/(1+exp(-p[2]/(kb*T)))
#  # gamma_up= p[0]*exp(-p[1]/(kb*T))*1/(1+exp((-p[1]+Ez)/(kb*T)))
#  # return p[0]*(1-exp(Ez/(kb*0.44)))/Ez
#  W=p[0]*Ez*Ez*(1/(1+exp((-p[1]+Ez)/(kb*T))))
#  return W

def Fermi(x,*p):
#         Vgs = np.array(eps)
    # alpha = p[1]
    # Vg0=p[2]
    e = 1.6e-19
    kb = 1.38e-23
    
#         dE = alpha * e * (Vgs-Vg0)*1e-3
    # dE = p[1] * e * (x-p[2])*1e-3


    fE = p[0]/(1+np.exp(-(p[1] * e * (x-p[2])*1e-3)/(kb*T)))
    return fE


popt=[]
pcov=[]
# popt,pcov = optimize.curve_fit(Fermi, x, y)
popt,pcov = curve_fit(Fermi,x,y,p0=[1,0.1,-770],maxfev=10000)

plot.plot(x,Fermi(x,*popt),label=r'alpha ='+str(np.round(popt[1],2))+', V0='+str(np.round(popt[2],2))+'mV')
plot.legend()
plot.savefig('FIT_ALPHA')


er=np.sqrt(np.diag(pcov))
erroralpha=er[1]



# popt[1]

plot.figure()
# dFermi=np.gradient(Fermi(x,*popt))
dFermi=np.gradient(y)
dFermi=phase=scipy.signal.savgol_filter(dFermi,11,5)
plot.plot(x,dFermi)
plt.savefig('Fermi_function')


