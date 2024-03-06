#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

#%%
import yrlib
import titanlib
import gridpp
import numpy as np
import netCDF4
from datetime import datetime

# Assemble data at this time point
frt = yrlib.util.date_to_unixtime(20231205,6)

# Variable
forecast_variable = "air_temperature_2m"

"""
        1. Import a single MEPS file
""" 

# Getting the data from yr
model = yrlib.model.Meps(include_ensemble=False)
dataset = model.get(frt, [forecast_variable],  
        start_time=frt, end_time= frt)


# The values of the forecast variable
meps_values = dataset.fields[forecast_variable].values
meps_values = meps_values[0,0,:,:]

# The corresponding grid points
meps_grid = dataset.gridpp_grid
meps_lats = meps_grid.get_lats()
meps_lons = meps_grid.get_lons()
meps_elevs = meps_grid.get_elevs()

#%%
"""
        2. Import a single NETATMO file
"""

if forecast_variable == "air_temperature_2m":
  short_variable = "ta"
elif forecast_variable == "relative_humidity_2m":
  short_variable = "uu"
netatmo_data = yrlib.netatmo.get([frt], short_variable) # ta = temperature, uu = humidity

# Unpack values, and convert to Kelvin
net_values = netatmo_data.values + 273.15 
net_values_0 = net_values[0,:] # Select only first slice (if there are several)
net_lons = netatmo_data.lons
net_lats = netatmo_data.lats
net_elevs = netatmo_data.elevs


# Removing na's
valid_indices = np.where((~np.isnan(net_values_0)))[0]

net_values_0 = net_values_0[valid_indices]
net_lons = net_lons[valid_indices]
net_lats = net_lats[valid_indices]
net_elevs = net_elevs[valid_indices]

#%%
"""
        3. Remove "invalid" netatmo observations with titanlib
"""

net_points = titanlib.Points(net_lats, net_lons, net_elevs)

# SCT settings
inner_radius = 50000
outer_radius = 100000
num_min = 5
num_max = 100
num_iterations = 1
num_min_prof = 20
dzmin = 100
dhmin = 10000
dz = 200
t2pos = np.full([len(net_values_0)], 4)
t2neg = np.full([len(net_values_0)], 4)
eps2 = np.full([len(net_values_0)], 0.5)

# Run the SCT, method to exclude bad observations
flags, sct, rep = titanlib.sct(net_points, net_values_0,
        num_min, num_max, inner_radius, outer_radius, num_iterations,
        num_min_prof, dzmin, dhmin , dz, t2pos, t2neg, eps2)


# The valid and invalid observations
index_valid_obs = np.where(flags == 0)[0]
index_invalid_obs = np.where(flags != 0)[0]
obsSize = index_valid_obs.shape[0]
obs_lats = net_lats[index_valid_obs]
obs_lons = net_lons[index_valid_obs]
obs_elevs = net_elevs[index_valid_obs]
obs_values = net_values_0[index_valid_obs]

# Date
yyyymmdd = yrlib.util.unixtime_to_date(frt)[0]
hh = yrlib.util.unixtime_to_date(frt)[1]
yyyy = int(yyyymmdd/10000)
mm = int((yyyymmdd-yyyy*10000)/100)
dd = yyyymmdd-yyyy*10000-100*mm
date = datetime(yyyy, mm, dd, hh, 0, 0)

# Write observations
with netCDF4.Dataset('observations.nc', 'w', format="NETCDF4") as file:
    nobs = file.createDimension('nobs', obsSize)
    ntime = file.createDimension('ntime', 6)
    nloc = file.createDimension('nloc', 3)
    ncol = file.createDimension('ncol', 2)
    time = file.createVariable('time',np.int32,('nobs', 'ntime'))
    loc = file.createVariable('loc',np.float64,('nobs', 'nloc'))
    cols = file.createVariable('cols',np.float64,('nobs', 'ncol'))
    cols.column_0 = "ObsVal"
    cols.column_1 = "ObsErr"
    for jo in range(0, obsSize):
        time[jo,:] = [date.year, date.month, date.day, date.hour, date.minute, date.second]
        loc[jo,:] = [obs_lons[jo], obs_lats[jo], obs_elevs[jo]]
        cols[jo,0] = obs_values[jo]
        cols[jo,1] = 0.1

# Write 5 observations for tests
with netCDF4.Dataset('observations_5.nc', 'w', format="NETCDF4") as file:
    nobs = file.createDimension('nobs', 5)
    ntime = file.createDimension('ntime', 6)
    nloc = file.createDimension('nloc', 3)
    ncol = file.createDimension('ncol', 2)
    time = file.createVariable('time',np.int32,('nobs', 'ntime'))
    loc = file.createVariable('loc',np.float64,('nobs', 'nloc'))
    cols = file.createVariable('cols',np.float64,('nobs', 'ncol'))
    cols.column_0 = "ObsVal"
    cols.column_1 = "ObsErr"
    for jo in range(0, 5):
        time[jo,:] = [date.year, date.month, date.day, date.hour, date.minute, date.second]
        loc[jo,:] = [obs_lons[jo], obs_lats[jo], obs_elevs[jo]]
        cols[jo,0] = obs_values[jo]
        cols[jo,1] = 0.1

#%%
"""
    4. Load finer grid (1km resolution)
"""

# Read fine grid from file
with netCDF4.Dataset('fine_grid.nc', 'r') as file:
    blats = file.variables['latitude'][:,:]
    blons = file.variables['longitude'][:,:]
    belevs = file.variables['altitude'][:,:]

# Downscaling
bgrid = gridpp.Grid(blats, blons, belevs)
background = gridpp.nearest(meps_grid, bgrid, meps_values)

# Write high-resolution background
with netCDF4.Dataset('background_hr.nc', 'w', format="NETCDF4") as file:
    nx = file.createDimension('nx', background.shape[1])
    ny = file.createDimension('ny', background.shape[0])
    nz = file.createDimension('nz_' + forecast_variable, 1)
    lats = file.createVariable('lats',np.float64,('ny','nx'))
    lons = file.createVariable('lons',np.float64,('ny','nx'))
    elevs = file.createVariable('elevs',np.float64,('ny','nx'))
    values = file.createVariable(forecast_variable,np.float64,('nz_' + forecast_variable,'ny','nx'))
    lats[:,:] = blats[:,:]
    lons[:,:] = blons[:,:]
    elevs[:,:] = belevs[:,:]
    values[0,:,:] = background[:,:]
