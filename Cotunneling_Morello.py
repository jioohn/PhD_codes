# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 18:31:51 2021

@author: AA255540
"""

##Morello model spin flip cotunneling

T=0.44

Ez=1.98*bohr*ydata

const=1

delta=1e-3


plt.figure()
delta=+400e-6


Gammares= 1/50e-6
# W0=kb*T*Gammares
W0=1

W=W0*Ez/np.sinh(Ez/(kb*T))

# tau=1/(const*(delta-Ez-220e-6)*(Ez/(1+exp(-Ez/(kb*T)))))



# W=const/(delta-Ez-220e-6)*(Ez/(1+exp(Ez/(kb*T))))
plt.figure()

plt.plot(Ez,1/W,label='cotunneling Delta= '+str(delta)+'eV')
plt.xlabel('Ez(eV)')
plt.ylabel('tau(A.U.)')
plt.legend()
plt.savefig('cotunnelingsimulation_gammavsB')




##### tunneling vs detuning


Bfield=0.1
T=0.44

Ez=1.98*bohr*Bfield
Ez=1.98*bohr*ydata

epsilon=220e-6


const=1

delta=450e-6
epsilon=220e-6

# epsilon=np.linspace(-1e-3,1e-3,201)

Gamma0=30e-6
W0=2*kb*T/h*(Gamma0/((delta+epsilon+Ez)*(delta+epsilon+Ez)))
deltaepsilon= np.sqrt(3*Gamma0*10e3/W0)

W=W0*Ez/np.sinh(Ez/(kb*T))/(1+(epsilon/deltaepsilon))


# tau=1/(const*(delta-Ez-220e-6)*(Ez/(1+exp(-Ez/(kb*T)))))

plt.figure()
delta=400e-6

# tau=const*(delta-Ez-220e-6)*(Ez/(1+exp(-Ez/(kb*T))))
plt.plot(epsilon,W,label='cotunneling Delta= '+str(delta)+'eV, B='+str(Bfield)+'T')
plt.xlabel('epsilon(eV)')
plt.ylabel('tau(A.U.)')
plt.legend()
plt.savefig('cotunnelingsimulation_tauvsepsilon')







############################
##fit cotunneling Morello

#####
# remove point at 0 Bfield

xaxis=xaxis[1:]
tempty=tempty[1:]
terr=terr[1:]
####



ydata=xaxis
T=0.44

Ez=1.98*bohr*ydata




plt.figure()
fig, ax = plt.subplots()
ax.set_xlabel('B (T)',fontsize=fontSize)
ax.set_ylabel(r'$ \tau_T$ (s)',fontsize=fontSize)
# plt.title(titolo)
#plt.plot(ydata,temptyfitarray_nodelay,'b+',label='data')
plt.errorbar(xaxis,tempty,yerr=terr,fmt='o',label='data from fit')
plt.legend(loc='upper left')
# plt.plot(xaxis,1e-8/Ez*np.sinh(Ez/(kb*T))-150e-6)



T=0.44
def cotunneling_Morello(x,*p):
 Ez=1.98*bohr*x


 t=p[0]*np.sinh((Ez+p[2])/(kb*T))/(Ez+p[2])
 return t +p[1]

# xaxis=ydata

popt,pcov = curve_fit(cotunneling_Morello,xaxis,tempty,p0=[1e-8,-100e-6,100e-6],maxfev=30000) 
err=er=np.sqrt(np.diag(pcov))

plt.plot(xaxis,cotunneling(xaxis,*popt),label='cotunneling formula')

plt.legend()
plt.savefig('fit_cotunnelingformula')

# Gamma_down=popt[3]/popt[0]

print('DeltaE='+str(popt[2]*1e6)+' microeV')
print('Gamma_lead2/pi*h='+str(popt[1])+'eV')
# 




# #######
# plt.figure()
# fig, ax = plt.subplots()
# ax.set_xlabel('B (T)',fontsize=fontSize)
# ax.set_ylabel(r'$ \tau_T$ (s)',fontsize=fontSize)
# # plt.title(titolo)
# #plt.plot(ydata,temptyfitarray_nodelay,'b+',label='data')
# plt.errorbar(xaxis,tempty,yerr=terr,fmt='o',label='data from fit')
# plt.legend(loc='upper left')
# # plt.plot(xaxis,1e-8/Ez*np.sinh(Ez/(kb*T))-150e-6)



# T=0.44
# def squared(x,*p):
#  Ez=1.98*bohr*x


#  t=p[0]*(1/Ez-p[1])**2+p[2]
#  return t 

# # xaxis=ydata

# popt,pcov = curve_fit(squared,xaxis,tempty,p0=[1,1,100e-6],maxfev=30000) 
# err=er=np.sqrt(np.diag(pcov))

# plt.plot(xaxis,squared(xaxis,*popt),label='cotunneling formula')

# plt.legend()
# plt.savefig('fit_cotunnelingformula')

# Gamma_down=popt[3]/popt[0]

# print('DeltaE='+str(popt[2]*1e6)+' microeV')
# print('Gamma_lead2/pi*h='+str(popt[1])+'eV')