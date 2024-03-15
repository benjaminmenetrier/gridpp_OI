#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

#%%
import numpy as np
import netCDF4

# Number of dirac points
ndx = 6
ndy = 6

# Variable
forecast_variable = "air_temperature_2m"

# Load background
#with netCDF4.Dataset('background_lr.nc', 'r') as file:
#    blons = file.variables['lon'][:,:]
#    blats = file.variables['lat'][:,:]
with netCDF4.Dataset('20231206_background_meps.nc', 'r') as file:
    blons = file.variables['lons'][:,:]
    blats = file.variables['lats'][:,:]

nx = blats.shape[1]
ny = blats.shape[0]
dx = nx/ndx
dy = ny/ndy

for idx in range(0,ndx):
  ix = int((0.5+idx)*dx)
  for idy in range(0,ndy):
    iy = int((0.5+idy)*dy)
    print("          {")
    print("            \"longitude\": \"" + str(blons[iy,ix]) + "\",")
    print("            \"latitude\": \"" + str(blats[iy,ix]) + "\",")
    print("            \"level\": \"1\",")
    print("            \"variable\": \"" + forecast_variable + "\"")
    print("          },")

