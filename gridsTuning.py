#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

#%%
import numpy as np
import netCDF4

# Config values
yamlLon = 17.50088882446289
yamlLat = 63.300415039062486

# Read fine grid from file
#with netCDF4.Dataset('fine_grid.nc', 'r') as file:
#    blats = file.variables['latitude'][:,:]
#    blons = file.variables['longitude'][:,:]
with netCDF4.Dataset('background_meps.nc', 'r') as file:
    blats = file.variables['lats'][:,:]
    blons = file.variables['lons'][:,:]

# Read ATLAS grid
#with netCDF4.Dataset('background_lr.nc', 'r') as file:
#    alats = file.variables['lat'][:,:]
#    alons = file.variables['lon'][:,:]
with netCDF4.Dataset('background_meps_test.nc', 'r') as file:
    alats = file.variables['lat'][:,:]
    alons = file.variables['lon'][:,:]

nx = alons.shape[0]
ny = alons.shape[1]
print("Center: ")
i = int(nx/2)
j = int(ny/2)
print(str(blons[i,j]) + " / " + str(blats[i,j]))
print(str(alons[i,j]) + " / " + str(alats[i,j]))
print("-----------------------")
print("Corner (0,0): ")
i = 0
j = 0
print(str(blons[i,j]) + " / " + str(blats[i,j]))
print(str(alons[i,j]) + " / " + str(alats[i,j]))
print("-----------------------")
print("Corner (nx-1,0): ")
i = nx-1
j = 0
print(str(blons[i,j]) + " / " + str(blats[i,j]))
print(str(alons[i,j]) + " / " + str(alats[i,j]))
print("-----------------------")
print("Corner (0,ny-1): ")
i = 0
j = ny-1
print(str(blons[i,j]) + " / " + str(blats[i,j]))
print(str(alons[i,j]) + " / " + str(alats[i,j]))
print("-----------------------")
print("Corner (nx-1,ny-1): ")
i = nx-1
j = ny-1
print(str(blons[i,j]) + " / " + str(blats[i,j]))
print(str(alons[i,j]) + " / " + str(alats[i,j]))
print("-----------------------")
norm = np.sqrt(np.average(np.square(blons-alons)+np.square(blats-alats)))
print("Norm: " + str(norm))
print("-----------------------")
print("New center:")
i = int(nx/2)
j = int(ny/2)
newLon = yamlLon+(blons[i,j]-alons[i,j])
newLat = yamlLat+(blats[i,j]-alats[i,j])
print("yamlLon = " + str(newLon))
print("yamlLat = " + str(newLat))
print("      \"lonlat(centre)\": [\"" + str(newLon) + "\",\"" + str(newLat) + "\"],")
