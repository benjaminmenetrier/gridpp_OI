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
#    blat = file.variables['latitude'][:,:]
#    blon = file.variables['longitude'][:,:]
with netCDF4.Dataset('background_meps.nc', 'r') as file:
    blat = file.variables['lat'][:,:]
    blon = file.variables['lon'][:,:]

# Read ATLAS grid
#with netCDF4.Dataset('background_lr.nc', 'r') as file:
#    alat = file.variables['lat'][:,:]
#    alon = file.variables['lon'][:,:]
with netCDF4.Dataset('background_meps_test.nc', 'r') as file:
    alat = file.variables['lat'][:,:]
    alon = file.variables['lon'][:,:]

nx = alon.shape[0]
ny = alon.shape[1]
print("Center: ")
i = int(nx/2)
j = int(ny/2)
print(str(blon[i,j]) + " / " + str(blat[i,j]))
print(str(alon[i,j]) + " / " + str(alat[i,j]))
print("-----------------------")
print("Corner (0,0): ")
i = 0
j = 0
print(str(blon[i,j]) + " / " + str(blat[i,j]))
print(str(alon[i,j]) + " / " + str(alat[i,j]))
print("-----------------------")
print("Corner (nx-1,0): ")
i = nx-1
j = 0
print(str(blon[i,j]) + " / " + str(blat[i,j]))
print(str(alon[i,j]) + " / " + str(alat[i,j]))
print("-----------------------")
print("Corner (0,ny-1): ")
i = 0
j = ny-1
print(str(blon[i,j]) + " / " + str(blat[i,j]))
print(str(alon[i,j]) + " / " + str(alat[i,j]))
print("-----------------------")
print("Corner (nx-1,ny-1): ")
i = nx-1
j = ny-1
print(str(blon[i,j]) + " / " + str(blat[i,j]))
print(str(alon[i,j]) + " / " + str(alat[i,j]))
print("-----------------------")
norm = np.sqrt(np.average(np.square(blon-alon)+np.square(blat-alat)))
print("Norm: " + str(norm))
print("-----------------------")
print("New center:")
i = int(nx/2)
j = int(ny/2)
newLon = yamlLon+(blon[i,j]-alon[i,j])
newLat = yamlLat+(blat[i,j]-alat[i,j])
print("yamlLon = " + str(newLon))
print("yamlLat = " + str(newLat))
print("      \"lonlat(centre)\": [\"" + str(newLon) + "\",\"" + str(newLat) + "\"],")
