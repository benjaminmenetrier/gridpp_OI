#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 Get and prepare raw data:
   1. Load MEPS background
   2. Read finer grid (1km resolution), downscale and write MEPS background
   3. Load and write observations 
"""

import argparse
import yrlib
import gridpp
import numpy as np
import netCDF4
from datetime import datetime

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("fineGrid", help="Fine grid file")
parser.add_argument("background", help="Background file")
parser.add_argument("observations", help="Observations file")
args = parser.parse_args()

"""
   Parameters
"""

# Variable
forecast_variable = "air_temperature_2m"

# Assemble data at this time point
frt = yrlib.util.date_to_unixtime(20231205,6)

"""
   1. Load MEPS background
"""

# Getting the data from yr
model = yrlib.model.Meps(include_ensemble=False)
dataset = model.get(frt, [forecast_variable], start_time=frt, end_time= frt)
meps_values = dataset.fields[forecast_variable].values[0,0,:,:]

# The corresponding grid points
meps_grid = dataset.gridpp_grid
meps_lat = meps_grid.get_lats()
meps_lon = meps_grid.get_lons()
meps_oro = meps_grid.get_elevs()

"""
   2. Load finer grid (1km resolution), downscale and write MEPS background
"""

# Read fine grid from file
with netCDF4.Dataset(args.fineGrid, 'r') as file:
  blat = file.variables['latitude'][:,:]
  blon = file.variables['longitude'][:,:]
  boro = file.variables['altitude'][:,:]

# Downscaling
bgrid = gridpp.Grid(blat, blon, boro)
background = gridpp.bilinear(meps_grid, bgrid, meps_values)

# Write background
with netCDF4.Dataset(args.background, 'w', format="NETCDF4") as file:
  nx = file.createDimension('nx', background.shape[1])
  ny = file.createDimension('ny', background.shape[0])
  nz_oro = file.createDimension('nz_oro', 1)
  nz_var = file.createDimension('nz_' + forecast_variable, 1)
  lat = file.createVariable('lat',np.float64,('ny','nx'))
  lon = file.createVariable('lon',np.float64,('ny','nx'))
  oro = file.createVariable('oro',np.float64,('nz_oro','ny','nx'))
  var = file.createVariable(forecast_variable,np.float64,('nz_' + forecast_variable,'ny','nx'))
  lat[:,:] = blat[:,:]
  lon[:,:] = blon[:,:]
  oro[0,:,:] = boro[:,:]
  var[0,:,:] = background[:,:]

"""
   3. Load and write observations 
"""

# Get NETATMO data
if forecast_variable == "air_temperature_2m":
  short_variable = "ta"
elif forecast_variable == "relative_humidity_2m":
  short_variable = "uu"
netatmo_data = yrlib.netatmo.get([frt], short_variable)

# Unpack values, and convert to Kelvin
net_values = netatmo_data.values[0,:] + 273.15 
net_lon = netatmo_data.lons
net_lat = netatmo_data.lats
net_oro = netatmo_data.elevs

# Removing na's
valid_indices = np.where((~np.isnan(net_values)))[0]
obs_values = net_values[valid_indices]
obs_lon = net_lon[valid_indices]
obs_lat = net_lat[valid_indices]
obs_oro = net_oro[valid_indices]
nobsAll = len(obs_lon)

# Define date components
yyyymmdd = yrlib.util.unixtime_to_date(frt)[0]
hh = yrlib.util.unixtime_to_date(frt)[1]
yyyy = int(yyyymmdd/10000)
mm = int((yyyymmdd-yyyy*10000)/100)
dd = yyyymmdd-yyyy*10000-100*mm
date = datetime(yyyy, mm, dd, hh, 0, 0)

# Write observations
with netCDF4.Dataset(args.observations, 'w', format="NETCDF4") as file:
  nobs = file.createDimension('nobs', nobsAll)
  ntime = file.createDimension('ntime', 6)
  nloc = file.createDimension('nloc', 3)
  ncol = file.createDimension('ncol', 1)
  time = file.createVariable('time',np.int32,('nobs', 'ntime'))
  loc = file.createVariable('loc',np.float64,('nobs', 'nloc'))
  cols = file.createVariable('cols',np.float64,('nobs', 'ncol'))
  cols.column_0 = "ObsVal"
  for jo in range(0, nobsAll):
    time[jo,:] = [date.year, date.month, date.day, date.hour, date.minute, date.second]
    loc[jo,:] = [obs_lon[jo], obs_lat[jo], obs_oro[jo]]
    cols[jo,0] = obs_values[jo]
