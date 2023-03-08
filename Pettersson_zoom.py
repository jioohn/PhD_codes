# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 17:19:31 2019

@author: AA255540
"""

zoom='S:\110-PHELIQS\110.05-LATEQS\110.05.01-QuantumSilicon\Tritonito\data\2019-04-12\mwfreqvsdetuning_ZOOM_5averages.txt'
with open(zoom,'r') as f:
     lines = f.readlines()
      #remove commas from datafile!!!
     z= np.array([line.split()[0] for line in lines],dtype='float64')
     
     dx=1.5
     dy=2.5
     eps=(1.5**2+2.5**2)**0.5
     x=np.linspace(0,eps,401)#mV
     y=-y
     y=y-y[0]
     yslope=(y[len(y)-1]-y[0])/len(y)