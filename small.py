#!/usr/bin/env python3

import os
import sys
from sys import exit
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import numpy as np
import cartopy.crs as ccrs

# Files
f_bkg = Dataset('background_small.nc', 'r', format='NETCDF4')
f_flat = Dataset('dirac_small_fastlam_flat.nc', 'r', format='NETCDF4')
f_bump = Dataset('dirac_small_bump.nc', 'r', format='NETCDF4')
f_fastlam = Dataset('dirac_small_fastlam.nc', 'r', format='NETCDF4')

# Coordinates
lon = f_bkg['lon'][:,:]
lat = f_bkg['lat'][:,:]
elevs = f_bkg['elevs'][0,:,:]

# Data
dirac_flat = f_flat['air_temperature_2m'][0,:,:]
dirac_bump = f_bump['air_temperature_2m'][0,:,:]
dirac_fastlam = f_fastlam['air_temperature_2m'][0,:,:]

# Levels
levels_elevs = np.linspace(0.0, 2250.0, 10)
levels_elevs_hr = np.linspace(0.0, 2250.0, 37)
levels_dirac = np.linspace(-1.0, 1.0, 21)

# Extent
extent = [5.3,9.1,60.2,62.1]

fig,ax = plt.subplots(ncols=2,nrows=2,figsize=(10,10),subplot_kw=dict(projection=ccrs.LambertConformal(central_longitude=7.2, central_latitude=61.15)))

ax[0][0].set_global()
ax[0][0].set_extent(extent, ccrs.PlateCarree())
ax[0][0].set_title("Orography")
ax[0][0].contourf(lon, lat, elevs, levels=levels_elevs_hr, cmap='terrain', transform=ccrs.PlateCarree())
ax[0][0].contour(lon, lat, elevs, levels=levels_elevs, linewidths=0.2, colors='k', transform=ccrs.PlateCarree())

ax[0][1].set_global()
ax[0][1].set_extent(extent, ccrs.PlateCarree())
ax[0][1].set_title("Increment without orography")
ax[0][1].contourf(lon, lat, elevs, levels=levels_elevs_hr, cmap='terrain', alpha=0.5, transform=ccrs.PlateCarree())
ax[0][1].contourf(lon, lat, dirac_flat, levels=levels_dirac, cmap='coolwarm', alpha=0.8, transform=ccrs.PlateCarree())

ax[1][0].set_global()
ax[1][0].set_extent(extent, ccrs.PlateCarree())
ax[1][0].set_title("Increment with orography - BUMP")
ax[1][0].contourf(lon, lat, elevs, levels=levels_elevs_hr, cmap='terrain', alpha=0.5, transform=ccrs.PlateCarree())
ax[1][0].contourf(lon, lat, dirac_bump, levels=levels_dirac, cmap='coolwarm', alpha=0.8, transform=ccrs.PlateCarree())

ax[1][1].set_global()
ax[1][1].set_extent(extent, ccrs.PlateCarree())
ax[1][1].set_title("Increment with orography - FastLAM")
ax[1][1].contourf(lon, lat, elevs, levels=levels_elevs_hr, cmap='terrain', alpha=0.5, transform=ccrs.PlateCarree())
ax[1][1].contourf(lon, lat, dirac_fastlam, levels=levels_dirac, cmap='coolwarm', alpha=0.8, transform=ccrs.PlateCarree())


plt.savefig("small.jpg", format="jpg", dpi=300)
plt.close()

os.system("mogrify -trim small.jpg")
