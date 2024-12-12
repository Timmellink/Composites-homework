# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 16:19:13 2024

@author: Timme
"""

# %% import statements
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
# %% import text file
widths = [6, 12, 14, 13]
df = pd.read_csv("3pb_data.txt", header=[0], sep='\s+')
#df2.columns = ["time", "disp", "force"]
# %% set variables

# %% plot data
df.plot(x="disp",y="force")
# %% calculate stiffness
# stiffness is 