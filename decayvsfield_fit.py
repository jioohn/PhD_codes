# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 17:52:09 2021

@author: AA255540
"""




#eliminate bad point for better fit @B=1.35


####

#account normalized probability
tempty=np.divide(temptyfitarray_nodelay,1000)
terr=np.divide(temptyerrorfitarray_nodelay,1000)
xaxis=ydata
####
# tempty=np.delete(tempty,27)
# terr=np.delete(terr,27)
# xaxis=np.delete(xaxis,27)




fig, ax = plt.subplots()
ax.set_xlabel('B (T)',fontsize=fontSize)
ax.set_ylabel(r'$ \tau_T$ (s)',fontsize=fontSize)
# plt.title(titolo)
#plt.plot(ydata,temptyfitarray_nodelay,'b+',label='data')
plt.errorbar(xaxis,tempty,yerr=terr,fmt='o',label='data from fit')
plt.legend(loc='upper left')




T=0.44
# T=0.44
def Pup_inv(x,*p):
 Ez=1.98*bohr*x
 E0=np.zeros(len(x))
 # gamma_down=p[0]*exp(-p[1]/(kb*T))*1/(1+exp(-p[1]/(kb*T)))
 # gamma_up= p[0]*exp(-p[1]/(kb*T))*1/(1+exp((-p[1]+Ez)/(kb*T)))
 gamma_down=exp(-p[0]*p[3])*1/(1+exp(-p[0]/(kb*T)))
 gamma_up= exp(-p[0]*p[4])*1/(1+exp((-p[0]+Ez)/(kb*T)))
 # return p[0]*(1-exp(Ez/(kb*0.44)))/Ez
 return p[2]*(gamma_up+gamma_down)/gamma_up +p[1]
# offset=-0.07
# tcoupling=5020e6*h+
 


parameter_bounds=([0,-30e-6,100e-6,-10e4,-10e4],[500e-6,30e-6,1e-3,10e4,10e4])
popt,pcov = curve_fit(Pup_inv,xaxis,tempty,p0=[100e-6,0,200e-6,1,1],bounds=parameter_bounds,maxfev=30000) 
err=er=np.sqrt(np.diag(pcov))

plt.plot(xaxis,Pup_inv(xaxis,*popt),label='gamma_down/gamma_up formula')

plt.legend()
plt.savefig('fit_doubleexp')

# Gamma_down=popt[3]/popt[0]

print('detuning='+str(popt[0])+'eV')
print('offset='+str(popt[1])+'s')
print('tau_0='+str(popt[2])+'s')
print('beta_down='+str(popt[3])+'1/eV')
print('beta_up='+str(popt[4])+'1/eV')
# print('Gamma_down='+str(popt[3])+'kHz')







##simulations Gamma vs B
kb=8.617333e-5#eV/K
beta=0.48#eV-1

detuning=44e-6


# fig, ax = plt.subplots()
gamma_down=exp(-detuning*beta)/(1+exp(-detuning/(kb*T)))
gamma_d=[]
for i in range(0,31):
    gamma_d.append(gamma_down)


gamma_up= exp(-detuning*beta)/(1+exp((-detuning+1.98*bohr*xaxis)/(kb*T)))
# plt.plot(xaxis,1/gamma_up,label='gamma_up')
# plt.plot(xaxis,1/gamma_d,label='gamma_down'ng
plt.plot(xaxis,(gamma_up+gamma_down)/(gamma_up),label='gamma_up')

plt.legend()








################################
fig, ax = plt.subplots()
ax.set_xlabel('B (T)',fontsize=fontSize)
ax.set_ylabel(r'$ \tau_T$ (s)',fontsize=fontSize)
# plt.title(titolo)
#plt.plot(ydata,temptyfitarray_nodelay,'b+',label='data')
plt.errorbar(xaxis,tempty,yerr=terr,fmt='o',label='data from fit')
plt.legend(loc='upper left')


T=0.44
def tauvsfield_4(x,*p):
 Ez=1.98*bohr*x
 # gamma_down=p[0]*exp(-p[2]/(kb*T))*1/(1+exp(-p[2]/(kb*T)))
 # gamma_up= p[0]*exp(-p[1]/(kb*T))*1/(1+exp((-p[1]+Ez)/(kb*T)))
 # return p[0]*(1-exp(Ez/(kb*0.44)))/Ez
 W=p[0]*Ez*Ez*(1/(1+exp((-p[1]+Ez)/(kb*T))))
 return W
# offset=-0.07
# tcoupling=500e6*h+
 
# T=0.355
xaxis=ydata

popt,pcov = curve_fit(tauvsfield_4,xaxis,tempty,p0=[1,100e-6],maxfev=30000) 
err=er=np.sqrt(np.diag(pcov))

plt.plot(xaxis,tauvsfield_4(xaxis,*popt),label='cotunneling formula')

plt.legend()
plt.savefig('fit_cotunnelingformula')

# Gamma_down=popt[3]/popt[0]

print('detuning='+str(popt[1]*1e6)+' microeV')
print('tau0='+str(popt[0])+'ms')





################################
xaxis=ydata
tempty=np.divide(temptyfitarray_nodelay,1000)
terr=np.divide(temptyerrorfitarray_nodelay,1000)

xaxis=xaxis[1:]
tempty=tempty[1:]
terr=terr[1:]
####



fig, ax = plt.subplots()
ax.set_xlabel('B (T)',fontsize=fontSize)
ax.set_ylabel(r'$ \tau_T$ (s)',fontsize=fontSize)
# plt.title(titolo)
plt.plot(xaxis,tempty,'b+',label='data')
# plt.errorbar(xaxis,tempty,yerr=terr,fmt='o',label='data from fit')
plt.legend(loc='upper left')


T=0.44
def cotunneling(x,*p):
 Ez=1.98*bohr*x
 # gamma_down=p[0]*exp(-p[2]/(kb*T))*1/(1+exp(-p[2]/(kb*T)))
 # gamma_up= p[0]*exp(-p[1]/(kb*T))*1/(1+exp((-p[1]+Ez)/(kb*T)))
 # return p[0]*(1-exp(Ez/(kb*0.44)))/Ez
 
 ###this formula works,  but  I am cheating with sign of Fermi function maybe
 W=p[0]/(p[1]+Ez)**2*(Ez/(1+exp(-Ez/(kb*T))))
 # W=p[0]/(p[1]+Ez)**2*(Ez/(1+exp(+Ez/(kb*T))))
 return 1/W+p[2]
# offset=-0.07
# tcoupling=500e6*h+
 
# T=0.355
# xaxis=ydata

popt,pcov = curve_fit(cotunneling,xaxis,tempty,p0=[1,100e-6,100e-6],maxfev=30000) 
err=er=np.sqrt(np.diag(pcov))

plt.plot(xaxis,cotunneling(xaxis,*popt),label='cotunneling formula')

plt.legend()
plt.savefig('fit_cotunnelingformula')

# Gamma_down=popt[3]/popt[0]

print('DeltaE='+str(popt[1]*1e6)+' microeV')
# print('Gamma_lead2/pi*h='+str(popt[1])+'eV')
# print('Gamma_lead2/pi='+str(popt[0])+'eV')







###simulation tau with cotunneling

T=0.44

Ez=1.98*bohr*ydata

const=1

delta=1e-3


plt.figure()
delta=+400e-6
W=const/((delta+Ez+220e-6)*(delta+Ez+220e-6))*(Ez/(1+exp(-Ez/(kb*T))))

W=const/((delta+Ez+220e-6)*(delta+Ez+220e-6))*(Ez/(1+exp(-Ez/(kb*T))))
# tau=1/(const*(delta-Ez-220e-6)*(Ez/(1+exp(-Ez/(kb*T)))))



# W=const/(delta-Ez-220e-6)*(Ez/(1+exp(Ez/(kb*T))))
plt.plot(Ez,W,label='cotunneling Delta= '+str(delta)+'eV')
plt.xlabel('Ez(eV)')
plt.ylabel('rate(A.U.)')
plt.legend()
plt.savefig('cotunnelingsimulation_tauvsB')



Bfield=0.1
T=0.44

Ez=1.98*bohr*Bfield

const=1

delta=1e-3
epsilon=np.linspace(-2e-3,2e-3,51)

tau=1/(const*(delta-Ez-epsilon)*(Ez/(1+exp(-Ez/(kb*T)))))

# tau=1/(const*(delta-Ez-220e-6)*(Ez/(1+exp(-Ez/(kb*T)))))

plt.figure()
delta=400e-6

# tau=const*(delta-Ez-220e-6)*(Ez/(1+exp(-Ez/(kb*T))))
plt.plot(epsilon,tau,label='cotunneling Delta= '+str(delta)+'eV, B='+str(Bfield)+'T')
plt.xlabel('epsilon(eV)')
plt.ylabel('tau(A.U.)')
plt.legend()
plt.savefig('cotunnelingsimulation_tauvsepsilon')





