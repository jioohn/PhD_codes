# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 13:16:07 2021

@author: AA255540
"""


import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp, linspace, random
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.special import factorial
from scipy.stats import poisson
import scipy

# # data = np.load(r'S:\132-PHELIQS\132.05-LATEQS\132.05.01-QuantumSilicon\Tritonito\2021\Feb2021_6G22_2_die103\exploration\G3=-772.429_g2=-997.286\time_traces\vG3f=0.0mV_vg2f=-0.0mV.npy')
# data = np.load(r'S:\132-PHELIQS\132.05-LATEQS\132.05.01-QuantumSilicon\Tritonito\2021\Feb2021_6G22_2_die103\exploration\G3=-773_g2=-995\time_traces\vG3f=-0.328mV_vg2f=0.051mV.npy')
# t=[]
# phi=[]
# for i  in range(0,len(data)):
#     t.append(data[i][0])
#     phi.append(data[i][1])
    
# plt.plot(t[0:1000],phi[0:1000])  




# data = np.load(r'S:\132-PHELIQS\132.05-LATEQS\132.05.01-QuantumSilicon\Tritonito\2021\Feb2021_6G22_2_die103\exploration\G3=-772.429_g2=-997.286\time_traces\vG3f=-0.6mV_vg2f=0.095mV.npy')
# t=[]
# phi=[]
# for i  in range(0,len(data)):
#     t.append(data[i][0])
#     phi.append(data[i][1])
    
# plt.plot(t,phi)  





# # vG3f=-0.312mV_vg2f=0.04mV

# data = np.load(r'S:/132-PHELIQS/132.05-LATEQS/132.05.01-QuantumSilicon/Tritonito/2021/Feb2021_6G22_2_die103/exploration/G3=-855_g2=-872/time_traces\vG3f=-0.312mV_vg2f=0.04mV.npy')
# t=[]
# phi=[]
# for i  in range(0,len(data)):
#     t.append(data[i][0])
#     phi.append(data[i][1])
# plt.plot(t,phi)    
# # plt.plot(t[4000:5000],phi[4000:5000])



# data = np.load(r'S:/132-PHELIQS/132.05-LATEQS/132.05.01-QuantumSilicon/Tritonito/2021/Feb2021_6G22_2_die103/exploration/G3=-855_g2=-872/time_traces\vG3f=-0.408mV_vg2f=0.044mV.npy')
# data = np.load(r'S:/132-PHELIQS/132.05-LATEQS/132.05.01-QuantumSilicon/Tritonito/2021/Feb2021_6G22_2_die103/exploration/G3=-855_g2=-872/time_traces\vG3f=-0.408mV_vg2f=0.044mV.npy')
# vg4f=-0.104mV_vg5f=0.02mV
# data=np.load(r'S:/132-PHELIQS/132.05-LATEQS/132.05.01-QuantumSilicon/Tritonito/2021/Feb2021_6G22_2_die103/exploration/g4=-736_g5=-1362/time_traces/vg4f=-0.128mV_vg5f=0.024mV.npy')
data=np.load(r'C:/Users/AA255540/Desktop/Tritonito2021_alldata/Feb2021_6G22_2_die103/exploration/g4=-736_g5=-1362/time_traces/vg4f=-0.128mV_vg5f=0.024mV.npy')


# data=np.load(r'S:/132-PHELIQS/132.05-LATEQS/132.05.01-QuantumSilicon/Tritonito/2021/Feb2021_6G22_2_die103/exploration/g4=-736_g5=-1362/time_traces/vg4f=-0.144mV_vg5f=0.027mV.npy')

# data=np.load(r'S:/132-PHELIQS/132.05-LATEQS/132.05.01-QuantumSilicon/Tritonito/2021/Feb2021_6G22_2_die103/exploration/g4=-736_g5=-1362/time_traces/vg4f=-0.064mV_vg5f=0.012mV.npy')

t=[]
phi=[]
for i  in range(0,len(data)):
    t.append(data[i][0])
    phi.append(data[i][1])
    
plt.figure()
plt.plot(t[0:500],phi[0:500])    
plt.xlabel('t(s)',fontSize=16)
plt.ylabel('$\phi_S(rad)$',fontSize=14)
plt.savefig('timetrace_5_vg4f=-0p128mV_vg5f=0p024mV')



######
Threshold=+0.07
T = t
PH = phi
tstep = T[1]-T[0]
Ttot = len(T)*tstep

# if Do_plot == True:
#     f = plt.figure()
#     f = plt.plot(T,PH)
Ts_up = [0]
Ts_down = [0]
#     t0 = time.time()
j = 0
for i,ph in enumerate(PH):
    if ph > Threshold:
        if i > 1 and PH[i-1] < Threshold:
            Ts_up.append(1)
            j = j+1
        elif i > 1 and PH[i-1] > Threshold:
            Ts_up[j] = Ts_up[j] + 1

j = 0
for i,ph in enumerate(PH):
    if ph < Threshold:
        if i > 1 and PH[i-1] < Threshold:
            Ts_down[j] = Ts_down[j] + 1  
        elif i > 1 and PH[i-1] > Threshold:
            Ts_down.append(1)
            j = j+1

P_up = sum(Ts_up)/len(PH)
P_down = sum(Ts_down)/len(PH)
T_up = np.mean(Ts_up)*tstep    # time spent in the two different states
T_down = np.mean(Ts_down)*tstep 


Ts_up=np.multiply(Ts_up,tstep)
Ts_down=np.multiply(Ts_down,tstep)



Ts_up=np.multiply(Ts_up,1e6)
Ts_down=np.multiply(Ts_down,1e6)

# def fit_Poissonian(x,*p):

#     '''poisson function, parameter lamb is the fit parameter'''
#       p[0]*p[1]^x*exp(p[1])/x
    
#     return p[0]*poisson.pmf(x, p[1])


def fit_Poissonian(x,*p):
    '''poisson function, parameter lamb is the fit parameter'''
    # f=(r*x)**y*exp(-r*x)/scipy.special.factorial(y)*Nevents
    return (p[0]*x)**p[1]*exp(-p[0]*x)/scipy.special.factorial(p[1])*p[2]





y,x=np.histogram(Ts_down,bins=50,range=(0,500))
# y=np.divide(y,np.sum(y))

x=x[0:-1]


popt,pcov = curve_fit(fit_Poissonian,x,y,p0=[1/100,1,834], maxfev=100000) # p0=[1e-3,PSmax,+0.014,phimin],bounds=parameter_bounds,


# plt.plot(y,fit_Poissonian(x,*popt),'+',label='fit')



plt.figure()
plt.hist(Ts_down,bins=50,range=(0,500),color='orange', label=r'$\tau_{in}$='+str(int(T_down*1e6))+'$\mu s$')

plt.xlabel(r'$\tau_{in}(\mu s)$',fontSize=14)
plt.ylabel('$counts$',fontSize=14)


mean_fit=np.sum(x*fit_Poissonian(x,*popt))/np.sum(fit_Poissonian(x,*popt))


plt.plot(x,fit_Poissonian(x,*popt),'black',label=r'fit $<\tau_{in}>$='+str(np.round(mean_fit))+'$\pm5$ $\mu s$, 1/r='+str(int(1/popt[0]))+'$\mu s$')
plt.legend()
plt.savefig('Tau_in_hist_fit_vg4f=-0.128mV_vg5f=0.024mV.png')








plt.figure()

y,x=np.histogram(Ts_up,bins=50,range=(0,500))
# y=np.divide(y,np.sum(y))

x=x[0:-1]


popt,pcov = curve_fit(fit_Poissonian,x,y,p0=[1/140,1,834], maxfev=100000) # p0=[1e-3,PSmax,+0.014,phimin],bounds=parameter_bounds,

err=np.diag(pcov)
# plt.plot(y,fit_Poissonian(x,*popt),'+',label='fit')



plt.plot()
plt.hist(Ts_up,bins=50,range=(0,500), label=r'$<\tau_{out}>$='+str(int(T_up*1e6))+'$\mu s$')

plt.xlabel(r'$\tau_{out}(\mu s)$',fontSize=14)
plt.ylabel('$counts$',fontSize=14)

mean_fit=np.sum(x*fit_Poissonian(x,*popt))/np.sum(fit_Poissonian(x,*popt))

# errmean=1/err[0]/np.sum(fit_Poissonian(x,*popt))


plt.plot(x,fit_Poissonian(x,*popt),'black',label=r'fit $<\tau_{out}>$='+str(int(mean_fit))+'$\pm31$ $\mu s$, 1/r='+str(int(1/popt[0]))+'$\mu s$')
plt.legend()
plt.savefig('Tau_out_hist_fit_vg4f=-0.128mV_vg5f=0.024mV.png')

# mean=np.mean(fit_Poissonian(x,*popt))




