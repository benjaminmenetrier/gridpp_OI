#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
  Compute errors with control observations
   1. Read errors
   2. Read background grid
   3. Plot errors
"""

import os
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import matplotlib as mp
import numpy as np
import cartopy.crs as ccrs
import netCDF4

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("errors", help="Errors")
parser.add_argument("background", help="Background")
args = parser.parse_args()

"""
   Parameters
"""

# Variable
forecast_variable = "air_temperature_2m"

"""
   1. Read errors
"""

# Read observations and background in observation space
with netCDF4.Dataset(args.errors, 'r') as file:
  lon = file.variables['lon'][:]
  lat = file.variables['lat'][:]
  bkgErr = file.variables['bkg'][:]
  oiErr = file.variables['oi'][:]
  varErr = file.variables['var'][:]

"""
   2. Read background grid
"""

# Read background grid
with netCDF4.Dataset(args.background, 'r') as file:
  blat = file.variables['lat'][:,:]
  blon = file.variables['lon'][:,:]

"""
   3. Plot errors
"""

# Compute error
diffErr = varErr-oiErr
diffErrNeg = np.count_nonzero(diffErr<0)/len(diffErr)

# Print statistics
print("Bkg:  " + str(np.mean(bkgErr)) + " +/- " + str(np.std(bkgErr)))
print("OI:   " + str(np.mean(oiErr)) + " +/- " + str(np.std(oiErr)))
print("Var:  " + str(np.mean(varErr)) + " +/- " + str(np.std(varErr)))
print("Diff: " + str(np.mean(diffErr)) + " +/- " + str(np.std(diffErr)) + " (" + str(100.0*diffErrNeg) + " %)")

# Domain bounds
lonMin = np.min(blon)
lonMax = np.max(blon)
latMin = np.min(blat)
latMax = np.max(blat)
central_longitude = 0.5*(lonMin+lonMax)
central_latitude = 0.5*(latMin+latMax)
extent = [lonMin,lonMax,latMin,latMax]

# Color map
cmAbs = mp.colormaps['jet']
cmRel = mp.colormaps['seismic']

# Bounds
vminAbs = 0.0
vmaxAbs = 0.0
vmaxAbs = max(vmaxAbs, np.max(bkgErr))
vmaxAbs = max(vmaxAbs, np.max(oiErr))
vmaxAbs = max(vmaxAbs, np.max(varErr))
vmaxRel = np.max(np.abs(diffErr))
vminRel = -vmaxRel

# Plot
fig,ax = plt.subplots(ncols=2,nrows=2,figsize=(10,10),subplot_kw=dict(projection=ccrs.LambertConformal(central_longitude=central_longitude, central_latitude=central_latitude)))

ax[0][0].set_global()
ax[0][0].coastlines(linewidth=0.2)
ax[0][0].set_extent(extent, ccrs.PlateCarree())
ax[0][0].set_title("Bkg error")
ax[0][0].scatter(lon, lat, s=10, c=bkgErr, cmap=cmAbs, vmin=vminAbs, vmax=vmaxAbs, transform=ccrs.PlateCarree())

ax[0][1].set_global()
ax[0][1].coastlines(linewidth=0.2)
ax[0][1].set_extent(extent, ccrs.PlateCarree())
ax[0][1].set_title("OI error")
ax[0][1].scatter(lon, lat, s=10, c=oiErr, cmap=cmAbs, vmin=vminAbs, vmax=vmaxAbs, transform=ccrs.PlateCarree())

ax[1][0].set_global()
ax[1][0].coastlines(linewidth=0.2)
ax[1][0].set_extent(extent, ccrs.PlateCarree())
ax[1][0].set_title("3dvar error")
ax[1][0].scatter(lon, lat, s=10, c=varErr, cmap=cmAbs, vmin=vminAbs, vmax=vmaxAbs, transform=ccrs.PlateCarree())

ax[1][1].set_global()
ax[1][1].coastlines(linewidth=0.2)
ax[1][1].set_extent(extent, ccrs.PlateCarree())
ax[1][1].set_title("3dvar-OI error")
ax[1][1].scatter(lon, lat, s=10, c=diffErr, cmap=cmRel, vmin=vminRel, vmax=vmaxRel, transform=ccrs.PlateCarree())

plt.show()
#plt.savefig("compareObservations.jpg", format="jpg", dpi=300)
#plt.close()

#os.system("mogrify -trim compareObservations.jpg")
