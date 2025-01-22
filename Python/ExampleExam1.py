# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 08:44:44 2025

@author: Timme
"""

#%% import statements
import numpy as np
import CompositeProperties as cp
import thermal_effects as te
import FailTest as ft

# %% set properties
layup = [0] + [45]*2+[90]*2+[-45]*2+[0]
h = 0.15e-3
Tcure = 200
Troom = 20
delta = Troom-Tcure
E1 = 100e9
E2 = 10e9
G12 = 5e9
nu12 = 0.3
alph1 = 0.2e-6
alph2 = 50e-6
S1t = 2100e6
S1c = 1400e6
S2t = 100e6
S2c = 120e6
S6 = 200e6