#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 21:18:07 2023

@author: StefaniaParigi
"""


kb=8.617333262145e-5 #eV/K
bohr=5.7883818012e-5 #eV/T
h=4.1357e-15  #eV*s

g=1.497

B=0.91

Ez=g*bohr*B
fl=Ez/h

print (fl*1e-9)