#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

#%%
import yrlib
import numpy as np
import netCDF4
from datetime import datetime

# Variable
forecast_variable = "air_temperature_2m"

# Obs source
pseudo_obs = False

# Date
frt = yrlib.util.date_to_unixtime(20231205,6)
yyyymmdd = yrlib.util.unixtime_to_date(frt)[0]
hh = yrlib.util.unixtime_to_date(frt)[1]
yyyy = int(yyyymmdd/10000)
mm = int((yyyymmdd-yyyy*10000)/100)
dd = yyyymmdd-yyyy*10000-100*mm
date = datetime(yyyy, mm, dd, hh, 0, 0)

# Load background
with netCDF4.Dataset('small_background.nc', 'r') as file:
  blat = file.variables['lat'][:,:]
  blon = file.variables['lon'][:,:]
  boro = file.variables['oro'][0,:,:]
  background = file.variables[forecast_variable][0,:,:]

if pseudo_obs:
  # Pseudo-observation location, value and date
  ix = [int(blon.shape[1]/2)+1, int(blon.shape[1]/2)-1]
  iy = [int(blon.shape[0]/2),int(blon.shape[0]/2)]
  delta = [1.0, -1.0]

  # Number of observations
  pseudoNobs = len(ix)

  # Write assimilated observations
  with netCDF4.Dataset('small_observations.nc', 'w', format="NETCDF4") as file:
    nobs = file.createDimension('nobs', pseudoNobs)
    ntime = file.createDimension('ntime', 6)
    nloc = file.createDimension('nloc', 3)
    ncol = file.createDimension('ncol', 2)
    time = file.createVariable('time',np.int32,('nobs', 'ntime'))
    loc = file.createVariable('loc',np.float64,('nobs', 'nloc'))
    cols = file.createVariable('cols',np.float64,('nobs', 'ncol'))
    cols.column_0 = "ObsVal"
    cols.column_1 = "ObsErr"
    for jobs in range(0, pseudoNobs):
      time[jobs,:] = [date.year, date.month, date.day, date.hour, date.minute, date.second]
      loc[jobs,:] = [blon[iy[jobs],ix[jobs]], blat[iy[jobs],ix[jobs]], boro[iy[jobs],ix[jobs]]]
      cols[jobs,:] = [background[iy[jobs],ix[jobs]]+delta[jobs], 1.0]
else:
  # Load observations
  with netCDF4.Dataset('observations.nc', 'r') as file:
    allLoc = file.variables['loc'][:,:]
    allCols = file.variables['cols'][:,:]

  # Get domain bounds
  lonMin = np.min(blon)
  lonMax = np.max(blon)
  latMin = np.min(blat)
  latMax = np.max(blat)

  # Keep observations inside the domain
  insideIobs = []
  for jobs in range(0, allLoc.shape[0]):
    if (lonMin <= allLoc[jobs,0]) and (allLoc[jobs,0] <= lonMax) and (latMin <= allLoc[jobs,1]) and (allLoc[jobs,1] <= latMax):
      insideIobs.append(jobs)

  # Number of observations
  insideNobs = len(insideIobs)

  # Write assimilated observations
  with netCDF4.Dataset('small_observations.nc', 'w', format="NETCDF4") as file:
    nobs = file.createDimension('nobs', insideNobs)
    ntime = file.createDimension('ntime', 6)
    nloc = file.createDimension('nloc', 3)
    ncol = file.createDimension('ncol', 2)
    time = file.createVariable('time',np.int32,('nobs', 'ntime'))
    loc = file.createVariable('loc',np.float64,('nobs', 'nloc'))
    cols = file.createVariable('cols',np.float64,('nobs', 'ncol'))
    cols.column_0 = "ObsVal"
    cols.column_1 = "ObsErr"
    for jobs in range(0, insideNobs):
      time[jobs,:] = [date.year, date.month, date.day, date.hour, date.minute, date.second]
      loc[jobs,:] = allLoc[insideIobs[jobs],:]
      cols[jobs,:] = allCols[insideIobs[jobs],:]
