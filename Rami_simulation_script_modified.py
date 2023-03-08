# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 18:37:54 2021

@author: AA255540
"""


import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

#filePath = 'data\\Simulations\\'
#fileName = 'cqNon0B'
#figC.savefig(filePath + fileName + '.png', transparent=True, bbox_inches='tight')

# Basis: S(0,2)
#        S(1,1)=|UpDown - DownUp>/sqrt(2), 
#        T+(1,1)=|UpUp>, 
#        T0(1,1)=|UpDown + DownUp>/sqrt(2)
#        T-(1,1)=|DownDown>

e = 1.60217662e-19
h = 4.135667662e-15 # in eV.s
muB = 5.778381e-5 # in eV/T
kB = 8.61733e-5 # in eV/K
g1 = 1.2 
g2 = 2
t = 20e-6 # in eV      NB: 1GHz = 4ueV
t0 = 1e-6 # in eV
tm = 7e-6 # in eV
tp = 0e-6 # in eV
sgn = -1
#tp and tm are inverted everywhere, is the same in the end tso

B = 0.6


Teff = t / kB # in K

Teff=0.1

Clim = 5

verticalLine = 0 # in ueV

numPointsE = 401
epsilonRange = np.linspace(-1e-4, 1e-4, numPointsE)
#epsilonRange = np.linspace(1.2e-4, 1.22e-4, numPointsE)

def update(val):
    global g1, g2, t, t0, tm, tp, sgn, Teff
    B = sB.val
    g1 = sG1.val
    g2 = sG2.val
    t = sT.val * 1e-6
    t0 = sT0.val * 1e-6
    tm = sTm.val * 1e-6
    tp = sTp.val * 1e-6
    Teff = sTeff.val * 2 * t / kB
    # solveEnergies(sgn, g1, g2, t, t0, tm, tp, B)
    solveEnergies(sgn, g1, g2, t, 0, tm, tm, B)
    for i in range(5):
            exec("e%d.set_ydata(np.multiply(diagram[%d], yMultiplier))" % (i, i))
            exec("c%d.set_ydata(np.multiply(cq[%d], 1e-4))" % (i, i))
    
    (figE.axes[0]).set_ylim(bottom=diagram.min() * yMultiplier, top=diagram.max() * yMultiplier)
    (figC.axes[0]).set_ylim(bottom=cq.min() * cMultiplier, top=cq.max() * cMultiplier)
    figE.canvas.draw_idle()
    figC.canvas.draw_idle()
    updateB(1)
    updateAvC()
    updateSpectro()

    
def flip(event):
    global sgn
    sgn = -sgn
    update(1)
    
def update2Ds(event):
    updateCB()
    updateCBBoltz()

def solveEnergies(sgn, g1, g2, t, t0, tm, tp, B):
    global energies, cq, diagram, ct
    energies = []
    cq = []
    ct = []
    for epsilon in epsilonRange:
        H = np.matrix([[-sgn/2*epsilon,       t,         tp,                                   t0,                        tm],
                        [t,             sgn/2*epsilon,    0,                               1/2*(g1-g2)*muB*B,               0],
                        [tp,                 0,         sgn/2*epsilon-1/2*(g1+g2)*muB*B,          0,                        0],
                        [t0,          1/2*(g1-g2)*muB*B,        0,                           sgn/2*epsilon,                 0],
                        [tm,                  0,                0,                             0,            sgn/2*epsilon+1/2*(g1+g2)*muB*B]])
        eigvals = (la.eigvalsh(H))
        energies.append(eigvals)
    diagram = ((np.array(energies)).reshape(numPointsE, 5)).transpose()
    for i in range(5):
        cq.append(-np.diff(np.diff(diagram[i])/np.diff(epsilonRange))/np.diff(epsilonRange)[:-1])
        ct.append(np.diff(diagram[i])/np.diff(epsilonRange))
    cq = ((np.array(cq)).reshape(5, numPointsE-2))    
    ct = ((np.array(ct)).reshape(5, numPointsE-1))    
    
# assuming sgn=-1
# def solveEnergies(sgn, g1, g2, t, t0, tm, tp, B):
#     global energies, cq, diagram, ct
#     energies = []
#     cq = []
#     ct = []
#     for epsilon in epsilonRange:
#         H = np.matrix([[-sgn/2*epsilon,       t,         tp,                                   t0,                        tm],
#                         [t,             sgn/2*epsilon,    0,                               1/2*(g1-g2)*muB*B,               0],
#                         [tp,                 0,         sgn/2*epsilon+1/2*(g1+g2)*muB*B,          0,                        0],
#                         [t0,          1/2*(g1-g2)*muB*B,        0,                           sgn/2*epsilon,                 0],
#                         [tm,                  0,                0,                             0,            sgn/2*epsilon-1/2*(g1+g2)*muB*B]])
#         eigvals = (la.eigvalsh(H))
#         energies.append(eigvals)
#     diagram = ((np.array(energies)).reshape(numPointsE, 5)).transpose()
#     for i in range(5):
#         cq.append(-np.diff(np.diff(diagram[i])/np.diff(epsilonRange))/np.diff(epsilonRange)[:-1])
#         ct.append(np.diff(diagram[i])/np.diff(epsilonRange))
#     cq = ((np.array(cq)).reshape(5, numPointsE-2))    
#     ct = ((np.array(ct)).reshape(5, numPointsE-1))    
    


fontSize = 20
xMultiplier = 1e6
yMultiplier = 1e6

#####
# t0=0
#####


solveEnergies(sgn, g1, g2, t, t0, tm, tp, B)

figE = plt.figure(1)
figE.clf()
#axE = figE.subplots()
plt.subplots_adjust(bottom=0.55)

plt.grid()
e0, = plt.plot(epsilonRange * xMultiplier, diagram[0] * yMultiplier, color='red')
e1, = plt.plot(epsilonRange * xMultiplier, diagram[1] * yMultiplier, color='gray')
e2, = plt.plot(epsilonRange * xMultiplier, diagram[2] * yMultiplier, color='blue')
e3, = plt.plot(epsilonRange * xMultiplier, diagram[3] * yMultiplier, color='brown')
e4, = plt.plot(epsilonRange * xMultiplier, diagram[4] * yMultiplier, color='green')
plt.xlabel('$\epsilon$ ($\mu eV$)', fontsize=fontSize)
plt.ylabel('E($\mu eV$)', fontsize=fontSize)
plt.tick_params(axis='both', which='major', labelsize=fontSize)

#plt.axvline(x=verticalLine, color='black')

axB = plt.axes([0.25, 0.2, 0.65, 0.03])
sB = Slider(axB, 'B (T)', 0, 1, valinit=B, valfmt='%1.2f',valstep=0.01)
sB.label.set_size(fontSize)
sB.valtext.set_fontsize(fontSize)
sB.on_changed(update)

axg1 = plt.axes([0.25, 0.25, 0.65, 0.03])
sG1 = Slider(axg1, '$g_1$', 1, 3, valinit=g1, valfmt='%1.2f',valstep=0.01)
sG1.label.set_size(fontSize)
sG1.valtext.set_fontsize(fontSize)
sG1.on_changed(update)

axg2 = plt.axes([0.25, 0.3, 0.65, 0.03])
sG2 = Slider(axg2, '$g_2$', 1, 3, valinit=g2, valfmt='%1.2f',valstep=0.01)
sG2.label.set_size(fontSize)
sG2.valtext.set_fontsize(fontSize)
sG2.on_changed(update)

axT = plt.axes([0.25, 0.35, 0.65, 0.03])
sT = Slider(axT, 't ($\mu eV$)', 0, 200, valinit=t*1e6, valfmt='%1.f',valstep=1)
sT.label.set_size(fontSize)
sT.valtext.set_fontsize(fontSize)
sT.on_changed(update)

# axT0 = plt.axes([0.25, 0.27, 0.65, 0.03])
# sT0 = Slider(axT0, '$t_0$ ($\mu eV$)', 0, 100, valinit=t0*1e6, valfmt='%1.f',valstep=1)
# sT0.label.set_size(fontSize)
# sT0.valtext.set_fontsize(fontSize)
# sT0.on_changed(update)

axTm = plt.axes([0.25, 0.4, 0.65, 0.03])
sTm = Slider(axTm, '$t_{SO}$ ($\mu eV$)', 0, 100, valinit=tm*1e6, valfmt='%1.f',valstep=1)
sTm.label.set_size(fontSize)
sTm.valtext.set_fontsize(fontSize)
sTm.on_changed(update)

# axTp = plt.axes([0.25, 0.37, 0.65, 0.03])
# sTp = Slider(axTp, '$t_p$ ($\mu eV$)', 0, 100, valinit=tp*1e6, valfmt='%1.f',valstep=1)
# sTp.label.set_size(fontSize)
# sTp.valtext.set_fontsize(fontSize)
# sTp.on_changed(update)

# axTeff = plt.axes([0.25, 0.42, 0.65, 0.03])
# sTeff = Slider(axTeff, '$T_{eff}$ (%2t)', 0.01, 3, valinit=Teff*kB/2/t, valfmt='%1.2f',valstep=0.01)
# sTeff.label.set_size(fontSize)
# sTeff.valtext.set_fontsize(fontSize)
# sTeff.on_changed(update)

recalcax = plt.axes([0.25, 0.02, 0.2, 0.04])
recalcButton = Button(recalcax, 'Recalculate', hovercolor='0.975')
recalcButton.on_clicked(update2Ds)

signax = plt.axes([0.7, 0.02, 0.2, 0.04])
signButton = Button(signax, 'Flip', hovercolor='0.975')
signButton.on_clicked(flip)

plt.draw()

# figE.savefig(r'Z:\132.05.01-QuantumSilicon\Triton experiments\Merit figures\Energy diagram spin qubit Square B 600mT without delta g.pdf', transparent=True, bbox_inches='tight')

    











































cMultiplier = 1e-4

figC = plt.figure(2)
figC.clf()
axC = figC.subplots()
plt.grid()
c0, = plt.plot(epsilonRange[:-2] * xMultiplier, cq[0] * cMultiplier, color='red')
c1, = plt.plot(epsilonRange[:-2] * xMultiplier, cq[1] * cMultiplier, color='gray')
c2, = plt.plot(epsilonRange[:-2] * xMultiplier, cq[2] * cMultiplier, color='blue')
c3, = plt.plot(epsilonRange[:-2] * xMultiplier, cq[3] * cMultiplier, color='brown')
c4, = plt.plot(epsilonRange[:-2] * xMultiplier, cq[4] * cMultiplier, color='green')
plt.xlabel('$\epsilon$ ($\mu eV$)', fontsize=fontSize)
plt.ylabel('$C_Q$ (a.u.)', fontsize=fontSize)
#plt.ylim(-15, 15)

#plt.axvline(x=verticalLine, color='black')

plt.draw()

# =============================================================================
# EDSR
# =============================================================================

epsilon = 0

d10 = []
d20 = []
d30 = []
d40 = []
d21 = []
d31 = []
d41 = []
d32 = []
d42 = []
d43 = []

numPointsB = 151
BRange = np.linspace(0, 1.5, numPointsB)

def resetDij():
    global d10, d20, d30, d40, d21, d31, d41, d32, d42, d43 
    d10 = []
    d20 = []
    d30 = []
    d40 = []
    d21 = []
    d31 = []
    d41 = []
    d32 = []
    d42 = []
    d43 = []
        
def solveEDSR(epsilon):
    resetDij()
    for B in BRange:
        H = np.matrix([[-sgn/2*epsilon, t, tp, t0, tm],
                       [t, sgn/2*epsilon, 0, 1/2*(g1-g2)*muB*B, 0],
                       [tp, 0, sgn/2*epsilon-1/2*(g1+g2)*muB*B, 0, 0],
                       [t0, 1/2*(g1-g2)*muB*B, 0, sgn/2*epsilon, 0],
                       [tm, 0, 0, 0, sgn/2*epsilon+1/2*(g1+g2)*muB*B]])
        eigvals = (la.eigvalsh(H))
    #    eigvals.sort()
        d10.append(abs(eigvals[1] - eigvals[0]))
        d20.append(abs(eigvals[2] - eigvals[0]))
        d30.append(abs(eigvals[3] - eigvals[0]))
        d40.append(abs(eigvals[4] - eigvals[0]))
        d21.append(abs(eigvals[2] - eigvals[1]))
        d31.append(abs(eigvals[3] - eigvals[1]))
        d41.append(abs(eigvals[4] - eigvals[1]))
        d32.append(abs(eigvals[3] - eigvals[2]))
        d42.append(abs(eigvals[4] - eigvals[2]))
        d43.append(abs(eigvals[4] - eigvals[3]))

def updateB(val):
    epsilon = sEps.val * 1e-6
    solveEDSR(epsilon)
    l10.set_xdata(np.multiply(d10, 1e-9/h))
    l20.set_xdata(np.multiply(d20, 1e-9/h))
    #l30.set_xdata(np.multiply(d%d%d, 1e-9/h))
    #l40.set_xdata(np.multiply(d%d%d, 1e-9/h))
    l21.set_xdata(np.multiply(d21, 1e-9/h))
    l31.set_xdata(np.multiply(d31, 1e-9/h))
    #l41.set_xdata(np.multiply(d%d%d, 1e-9/h))
    #l32.set_xdata(np.multiply(d%d%d, 1e-9/h))
    #l42.set_xdata(np.multiply(d%d%d, 1e-9/h))
    #l43.set_xdata(np.multiply(d%d%d, 1e-9/h))
    figB.canvas.draw_idle()
    figE.canvas.draw_idle()



fontSize = 16
figB = plt.figure(3)
figB.clf()
ax = figB.subplots()
plt.subplots_adjust(bottom=0.25)


solveEDSR(epsilon)
plt.grid()
l10, = plt.plot(np.multiply(d10, 1e-9/h), BRange, color='red')
l20, = plt.plot(np.multiply(d20, 1e-9/h), BRange, color='blue')
#l30, = plt.plot(np.multiply(d30, 1e-9/h), BRange, color='green')
#l40, = plt.plot(np.multiply(d40, 1e-9/h), BRange, color='brown')
l21, = plt.plot(np.multiply(d21, 1e-9/h), BRange, color='black')
l31, = plt.plot(np.multiply(d31, 1e-9/h), BRange, color='gray')
#l41, = plt.plot(np.multiply(d41, 1e-9/h), BRange, color='pink')
#l32, = plt.plot(np.multiply(d32, 1e-9/h), BRange, color='purple')
#l42, = plt.plot(np.multiply(d42, 1e-9/h), BRange, color='yellow')
#l43, = plt.plot(np.multiply(d43, 1e-9/h), BRange, color='orange')

plt.xlabel('Frequency (GHz)', fontsize=fontSize)
plt.ylabel('B (T)', fontsize=fontSize)
plt.tick_params(axis='both', which='major', labelsize=fontSize)

axEps = plt.axes([0.25, 0.05, 0.65, 0.03])
sEps = Slider(axEps, 'Detuning ($\mu eV$)', -100, 100, valinit=epsilon, valfmt='%1.f',valstep=1)
sEps.label.set_size(fontSize)
sEps.valtext.set_fontsize(fontSize)
sEps.on_changed(updateB)

plt.draw()
plt.show()

# =============================================================================
# Interdot vs B
# =============================================================================

def solveCB():
    global cqb
    cqb = []
    for B in BRange:
        ground = []
        for epsilon in epsilonRange:
            H = np.matrix([[-sgn/2*epsilon, t, tp, t0, tm],
                           [t, sgn/2*epsilon, 0, 1/2*(g1-g2)*muB*B, 0],
                           [tp, 0, sgn/2*epsilon-1/2*(g1+g2)*muB*B, 0, 0],
                           [t0, 1/2*(g1-g2)*muB*B, 0, sgn/2*epsilon, 0],
                           [tm, 0, 0, 0, sgn/2*epsilon+1/2*(g1+g2)*muB*B]])
            eigvals = (la.eigvalsh(H))
            ground.append(eigvals[0])
        cqb.append(-np.diff(np.diff(ground)/np.diff(epsilonRange))/np.diff(epsilonRange)[:-1])
    cqb = ((np.array(cqb)).reshape(numPointsB, numPointsE-2))

solveCB()
figCB = plt.figure(4)
figCB.clf()

CvB = plt.pcolormesh(epsilonRange[:-2] * xMultiplier, BRange, cqb * cMultiplier, cmap='plasma')

plt.xlabel('$\epsilon$ ($\mu eV$)', fontsize=fontSize)
plt.ylabel('B (T)', fontsize=fontSize)
plt.tick_params(axis='both', which='major', labelsize=fontSize)
cb = plt.colorbar()
cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=fontSize)
cb.set_label(label='$C_Q$ (a.u.)', fontsize=fontSize)


plt.draw()
plt.show()

def updateCB():
    global CvB, cb, figCB
    solveCB()
    plt.figure(4)
    CvB = plt.pcolormesh(epsilonRange[:-2] * xMultiplier, BRange, cqb * cMultiplier, cmap='plasma')
    cb.on_mappable_changed(CvB)
    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=fontSize)
    figCB.canvas.draw_idle()
    
# =============================================================================
# Boltzmann
# =============================================================================

def solveBoltzmann():
    global diagram, probs, probabilities, cq, averagedCq, ct, wheightedCt, cTot
    probs = []
    beta = 1/(kB*Teff)
    for i in range(numPointsE):
        exps = np.exp(-beta * diagram[:,i])
        Z = exps.sum()
        probs.append(exps/Z)
    probabilities = ((np.array(probs)).reshape(numPointsE, 5)).transpose()
    averagedCq = (np.multiply(probabilities[:,:-2], cq)).sum(axis=0)
    dPdEps = np.divide(np.diff(probabilities), np.diff(epsilonRange))
    wheightedCt = (np.multiply(dPdEps[:,:-1], ct[:,:-1])).sum(axis=0)
    cTot = averagedCq + wheightedCt
    

solveBoltzmann()
figAvC = plt.figure(5)
figAvC.clf()
plt.grid()
avC, = plt.plot(epsilonRange[:-2] * xMultiplier, averagedCq * cMultiplier, color='red')
whC, = plt.plot(epsilonRange[:-2] * xMultiplier, wheightedCt * cMultiplier, color='blue')
totC, = plt.plot(epsilonRange[:-2] * xMultiplier, cTot * cMultiplier, color='purple')
plt.xlabel('$\epsilon$ ($\mu eV$)', fontsize=fontSize)
plt.ylabel('$C_Q$ (a.u.)', fontsize=fontSize)
#plt.axvline(x=verticalLine, color='black')
plt.draw()

def updateAvC():
    solveBoltzmann()
    newAv = averagedCq * cMultiplier
    newWh = wheightedCt * cMultiplier
    newTot = cTot * cMultiplier
    avC.set_ydata(newAv)
    whC.set_ydata(newWh)
    totC.set_ydata(newTot)
    plt.figure(5)
    plt.ylim(min(newAv.min(), newWh.min(), newTot.min()), max(newAv.max(), newWh.max(), newTot.max()))
    figAvC.canvas.draw_idle()
    
def solveCbBoltzmann():
    global cqBo, cTotBo
    cqBo = []
    cTotBo = []
    for B in BRange:
        solveEnergies(sgn, g1, g2, t, t0, tm, tp, B)
        solveBoltzmann()
        cqBo.append(averagedCq)
        cTotBo.append(cTot)
    cqBo = (np.array(cqBo)).reshape(numPointsB, numPointsE-2)
    cTotBo = (np.array(cTotBo)).reshape(numPointsB, numPointsE-2)
    
solveCbBoltzmann()
figCBBoltz = plt.figure(6)
figCBBoltz.clf()

CvBBoltz = plt.pcolormesh(epsilonRange[:-2] * xMultiplier, BRange, cqBo * cMultiplier, cmap='plasma')

plt.xlabel('$\epsilon$ ($\mu eV$)', fontsize=fontSize)
plt.ylabel('B (T)', fontsize=fontSize)
plt.tick_params(axis='both', which='major', labelsize=fontSize)
cbB = plt.colorbar()
cbB.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=fontSize)
cbB.set_label(label='$C_Q$ (a.u.)', fontsize=fontSize)


plt.draw()
plt.show()

figCTotBBoltz = plt.figure(7)
figCTotBBoltz.clf()

CTotvBBoltz = plt.pcolormesh(epsilonRange[:-2] * xMultiplier, BRange, cTotBo * cMultiplier, cmap='plasma')

plt.xlabel('$\epsilon$ ($\mu eV$)', fontsize=fontSize)
plt.ylabel('B (T)', fontsize=fontSize)
plt.tick_params(axis='both', which='major', labelsize=fontSize)
cbBTot = plt.colorbar()
cbBTot.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=fontSize)
cbBTot.set_label(label='$C_{Tot}$ (a.u.)', fontsize=fontSize)


plt.draw()
plt.show()

def updateCBBoltz():
    global CvBBoltz, cbB, figCBBoltz
    solveCbBoltzmann()
    plt.figure(6)
    CvBBoltz = plt.pcolormesh(epsilonRange[:-2] * xMultiplier, BRange, cqBo * cMultiplier, cmap='plasma')
    cbB.on_mappable_changed(CvBBoltz)
    cbB.ax.set_yticklabels(cbB.ax.get_yticklabels(), fontsize=fontSize)
    figCBBoltz.canvas.draw_idle()
    plt.figure(7)
    CTotvBBoltz = plt.pcolormesh(epsilonRange[:-2] * xMultiplier, BRange, cTotBo * cMultiplier, cmap='plasma')
    cbBTot.on_mappable_changed(CTotvBBoltz)
    cbBTot.ax.set_yticklabels(cbBTot.ax.get_yticklabels(), fontsize=fontSize)
    figCTotBBoltz.canvas.draw_idle()
    
# =============================================================================
# Spectroscopy
# =============================================================================

figS = plt.figure(8)
figS.clf()

plt.grid()
s10, = plt.plot(epsilonRange * xMultiplier, np.multiply(diagram[2] - diagram[0], 1e-9/h), color='red')

plt.xlabel('$\epsilon$ ($\mu eV$)', fontsize=fontSize)
plt.ylabel('Frequency (GHz)', fontsize=fontSize)
plt.tick_params(axis='both', which='major', labelsize=fontSize)

plt.draw()
plt.show()

def updateSpectro():
    global figS, diagram, s10
    newData = np.multiply(diagram[2] - diagram[0], 1e-9/h)
    s10.set_ydata(newData)
    plt.figure(7)
    plt.ylim(newData.min(), newData.max())
    figS.canvas.draw_idle()