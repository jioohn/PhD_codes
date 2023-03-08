# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 16:08:28 2019

@author: AA255540
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import exp, linspace, random
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from scipy.signal import savgol_filter
import math

#no FIT just smooth it an extract extract FWHM 
#
def lorentz(x,*p):
   return p[0] /((3.14159/2)*p[2])* 1/(1.0 + ((x  - p[1])/p[2])**2)

#treat 2d data
    #     z= np.array([line.split()[2] for line in lines[3:50]],dtype='float64')
    #     xdata=x
    #     ydata=y
    #     zdata=z
    #     avg=zdata
    ##    dx=x[0]-x[1]
    ##    dy = np.diff(y)/dx
kb=0.000086167#(eV/K)
T=0.44
beta=1/(kb*T)
alpha=0.25

def interdotfitT0(x,*t):
     return t[1]*(2*t[0])**2 / ( (2*t[0])**2+ (x*alpha)**2)**1.5 
def interdotfit_alphafit(x,*p):
    #p[0]-->t
    #p[1]-->alpha
    #p[2]--> E
    #p[3-->const=1/Z
#    p[4]-->T
#    return (2*t)**2 / ( (2*t)**2+ x**2)**1.5*exp(-E/(kb*T))*const
     return (2*p[0])**2 / ( (2*p[0])**2+ (p[1]*x)**2)**1.5 *exp( ( (2*p[0])**2+ (p[1]*x)**2)**0.5 *beta)*p[2]

def interdotfit_tempfit(x,*p):
    #p[0]-->t
    #p[1]-->alpha
    #p[2]--> E
    #p[3-->const=1/Z
#    p[4]-->T
#    return (2*t)**2 / ( (2*t)**2+ x**2)**1.5*exp(-E/(kb*T))*const
     return (2*p[0])**2 / ( (2*p[0])**2+ (p[1]*x)**2)**1.5 *exp( ( (2*p[0])**2+ (p[1]*x)**2)**0.5 *p[3])*p[2]
 
def interdotfit_alphafix(x,*p):
    #p[0]-->t
    #p[1]-->const
#    return (2*t)**2 / ( (2*t)**2+ x**2)**1.5*exp(-E/(kb*T))*const
     return (2*p[0])**2 / ( (2*p[0])**2+ (alpha*x)**2)**1.5 *exp( ( (2*p[0])**2+ (alpha*x)**2)**0.5 *beta)*p[1]   


h=4.135667*1e-15
Nmeas=8
x = []
y = []
z=[]
xdata=[]
ydata=[]
zdata=[]
xsimm=[]
name=['']*Nmeas
#filename=S:\\101-PHELIQS\\101.05-LATEQS\\101.05.01-QuantumSilicon\\Tritonito2\\resume_lowQ.txt
folder= r'S:\110-PHELIQS\110.05-LATEQS\110.05.01-QuantumSilicon\Tritonito\data\2019-04-05'
xpoints=401
ypoints=121
startdata=3
width=[]*Nmeas
erwidth=[]*Nmeas
fwhm=[]*Nmeas
maxphi=[]*Nmeas
ermax=[]*Nmeas
fitT0=[0]*xpoints
guessT=[0]*xpoints
fontSize=16
#name[0]='\TP4_T410mK_Tdependance_notaveraged.txt'#better one averaged
name[0]='\TP4_T490mK_Tdependance_averaged.txt'
name[1]='\TP4_T540mK_Tdependance_averaged.txt'
name[2]='\TP4_T580mK_Tdependance_averaged.txt'
name[3]='\TP4_T650mK_Tdependance_averaged.txt'
name[4]='\TP4_T700mK_Tdependance_averaged.txt'
name[5]='\TP4_T800mK_Tdependance_averaged.txt'
name[6]='\TP4_T1K_Tdependance_averaged.txt'
name[7]='\TP4_T1p5K_Tdependance_averaged.txt'
tcoupling=[]
errt=[]
Tax=[490,540,580,650,700,800,1000,1500]   

for i in range(0,Nmeas):

    
    with open(folder+name[i],'r') as f:
     lines = f.readlines()
      #remove commas from datafile!!!
     y= np.array([line.split()[0] for line in lines],dtype='float64')
     dx=1.5
     dy=2.5
     eps=(1.5**2+2.5**2)**0.5
     x=np.linspace(0,eps,401)#mV
     y=-y
     y=y-y[0]
     yslope=(y[len(y)-1]-y[0])/len(y)
#     t=5e9#GHz
#     t=t*h

     #the x=0 should be on the center of the peak
     
#     alpha=0.08
#     xsimm=xsimm*alpha #eV
     #insertingparam manual
#     t=3e9 #(GHz)
#     t=h*t#eV
#     alpha=0.15
#         
     T=Tax[i]
     for k in range(0,len(y)):
         y[k]=y[k]-yslope*k
         #guess of the curve
#         guessT[k] =(2*t)**2 / ( (2*t)**2+ (alpha*xsimm[k])**2)**1.5 *exp( -( (2*t)**2+ (alpha*xsimm[k])**2)**0.5 *beta)#better extract parameter from fit
    
#    ysmooth=savgol_filter(y, 51, 3)   
     ysmooth=y    
#    ysmooth=ysmooth*1000
    x0=np.argmax(ysmooth)
    deltax=x[x0]-x[len(x)-1]/2
    xsimm=linspace(-x[len(x)-1]/2,x[len(x)-1]/2, len(x))
    xsimm=xsimm-deltax#in meV
    xsimm=xsimm/1000#V   



    fig, ax = plt.subplots()
 
    ax.set_xlabel('$\epsilon$ (mV)',fontsize=fontSize)
    ax.set_ylabel('$\phi $(mrad)',fontsize=fontSize)
#    plt.plot(xsimm,guessT/max(guessT))
    plt.title('Fit at T='+str(T)+'mK')
    plt.plot(xsimm*1000,y*1000,'b+',label='data@'+str(T)+'mK')#mV on x  
    plt.show()
    
    

    alpha=0.23
    beta=1/(kb*T)
#   E=0.0001#100ueV
    t=10*1e9#GHz
    t=h*t
    const=0.1#const=alpha**2/Z#Z partition function
    offset=0
    def interdotfit_alphafix(x,*p):
    #p[0]-->t
    #p[1]-->const
#    return (2*t)**2 / ( (2*t)**2+ x**2)**1.5*exp(-E/(kb*T))*const
     return (2*p[0])**2 / ( (2*p[0])**2+ (alpha*x)**2)**1.5 *exp( -( (2*p[0])**2+ (alpha*x)**2)**0.5 *beta/2)*p[1] +p[2]
   
    def interdotfit_alphafit(x,*p):
    #p[0]-->t
    #p[1]-->alpha
    #p[2]--> const
#    return (2*t)**2 / ( (2*t)**2+ x**2)**1.5*exp(-E/(kb*T))*const
     return (2*p[0])**2 / ( (2*p[0])**2+ (p[1]*x)**2)**1.5 *exp( -( (2*p[0])**2+ (p[1]*x)**2)**0.5 *beta/2)*p[2] 
    
#    popt,pcov = curve_fit(interdotfit_alphafit,xsimm,ysmooth,p0=[t,alpha,const],maxfev=2000)
    
    popt,pcov = curve_fit(interdotfit_alphafix,xsimm,ysmooth,p0=[t,const,offset],maxfev=2000)
#    plt.plot(xsimm,fitT0/max(fitT0),'r',label='guess_renormalized')
    norm=max(interdotfit_alphafix(xsimm,*popt))
    plt.plot(xsimm*1000,interdotfit_alphafix(xsimm,*popt)*1000,'r',label='fit@'+str(T)+', alpha='+str(alpha))
    
#    norm=max(interdotfit_alphafit(xsimm,*popt))
#    plt.plot(xsimm,interdotfit_alphafit(xsimm,*popt)/norm,'r',label='fit 440mK, $)
    tfit=popt[0]/h*1e-9
    terrorfit=pcov[0]/h*1e-9
    er=np.sqrt(np.diag(pcov))
    errort=er[0]/h*1e-9
    errt.append(float(errort))
    
    tcoupling.append(tfit)#(GHZ)

    print('t='+str(tfit)+'+- '+str(errort)+' GHz, alpha='+str(alpha))
    plt.legend(loc='upper left')
    fig.savefig(folder+'TP4_'+str(T)+'mK_10averages_nosmooth_slopeadjust_V_rad_fit_fixedalpha_0p23')






##fit



#    xax=np.divide(xax,1000) #K
#    T=xax[i]
#    beta=kb*T
#    def interdotfit(x,*p):
#     return p[2]*(2*p[0])**2 / ( (2*p[0])**2+ (p[1]*x)**2)**1.5 *exp(-E/(beta))
#
#    popt,pcov = curve_fit(interdotfit,xsimm,ysmooth,p0=[t,const,alpha],maxfev=1000)
#    
#    #    def interdotfit(x,*p):
##     return (2*p[0])**2 / ( (2*p[0])**2+ (p[1]*x)**2)**1.5 *exp(-p[2]/(kb*p[4]))*p[3]
#    popt,pcov = curve_fit(interdotfit,xsimm,ysmooth,p=[t,alpha,E,const,T],maxfev=10000)
#
#
#    plt.plot(xsimm,ysmooth,'b+:',label='data')
#    plt.plot(xsimm,interdotfit(xsimm,*popt),'ro:',label='fit')
#    plt.plot(xsimm,fitT0)

    

    
    phimax=norm
    kmin=0
    kmax=0
    for k in range(0,len(ysmooth)):
        if(interdotfit_alphafix(xsimm,*popt)[k]>phimax/2):
            kmax=k
            if kmin==0:
                kmin=k
    maxphi.append(norm*1000)      
    width=(kmax-kmin)*eps/401 
    fwhm.append(width)    
    

#    erwidth.append(2.355*perr[2])

    

 
fig, ax = plt.subplots()   
plt.plot(Tax,fwhm,'b+')
plt.xlim(400,1550)
plt.xlabel('$T$ (mK)')
plt.ylabel('$FWHM$ (mV)')
#plt.errorbar(Tax, fwhm ,yerr=erwidth)
fig.savefig(folder+'FWHMvsT') 


fig, ax = plt.subplots()   
#plt.plot(Tax,tcoupling,'b+')
plt.xlabel('$T$ (mK)')
plt.ylabel('$t$ (GHz)')
ax.errorbar(Tax,tcoupling ,yerr=errt,fmt='o')
fig.savefig(folder+'tcouplingvsT') 



fig, ax = plt.subplots()   
plt.plot(Tax,maxphi,'b+')
plt.xlim(400,1550)
plt.xlabel('$T$ (mK)')
plt.ylabel('${\phi_{min}}$ (mrad)')
fig.savefig(folder+'max(phi)_fromfit_vsT') 
#plt.errorbar(Tax, width,yerr=maxer)






t=linspace(1,1e-3,100)#eV
epsilon=linspace(-3e-3,3e-3,101)
T=0.44
alpha=1
beta=1/(kb*T)
alphaprime=1#alpha1-alpha2
electron=1.6e-19
x=epsilon
#interdot simulation
for i in range(0,101):
    t=t[0]
    f=((alphaprime*electron)**2)*(2*t)**2 / ( 2*(t)**2+ alpha*x**2)**1.5*exp(-( (2*t)**2+ (alpha*x)**2)**0.5 *beta/2)


def extractFWHM(f,epsilon):
    kmin=0
    kmax=0
    phimax=max(f)
    epsrange=min(eps)-max(eps)
    for k in range(0,len(epsilon)):
            if(f[k]>phimax/2):
                kmax=k
                if kmin==0:
                    kmin=k
                       
    width=(kmax-kmin)*epsrange/len(epsilon) 
return width phimax