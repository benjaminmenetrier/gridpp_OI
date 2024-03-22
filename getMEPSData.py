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

# Variable
forecast_variable = "air_temperature_2m"

# List of cycles
yyyymmdd_list = [20231126,20231127,20231128,20231129,20231130,20231201,20231202,20231204,20231205,20231206]

for yyyymmdd in yyyymmdd_list:
  # Assemble data at this time point
  frt = yrlib.util.date_to_unixtime(yyyymmdd,6)

  # Getting the data from yr
  model = yrlib.model.Meps(include_ensemble=True)
  dataset = model.get(frt, [forecast_variable], start_time=frt, end_time= frt)
  meps_values = dataset.fields[forecast_variable].values[0,0:15,:,:]

  # The corresponding grid points
  meps_grid = dataset.gridpp_grid
  meps_lat = meps_grid.get_lats()
  meps_lon = meps_grid.get_lons()
  meps_oro = meps_grid.get_elevs()

  # Write meps background (ensemble mean)
  meps_avg = np.average(meps_values, axis=0)
  with netCDF4.Dataset(str(yyyymmdd) + '_background_meps.nc', 'w', format="NETCDF4") as file:
    nx = file.createDimension('nx', meps_avg.shape[1])
    ny = file.createDimension('ny', meps_avg.shape[0])
    nz = file.createDimension('nz_' + forecast_variable, 1)
    lat = file.createVariable('lat',np.float64,('ny','nx'))
    lon = file.createVariable('lon',np.float64,('ny','nx'))
    oro = file.createVariable('oro',np.float64,('ny','nx'))
    values = file.createVariable(forecast_variable,np.float64,('nz_' + forecast_variable,'ny','nx'))
    lat[:,:] = meps_lat[:,:]
    lon[:,:] = meps_lon[:,:]
    oro[:,:] = meps_oro[:,:]
    values[0,:,:] = meps_avg[:,:]

  for ie in range(0, meps_values.shape[0]):
    # Write meps members
    with netCDF4.Dataset(str(yyyymmdd) + '_member_meps_' + str(ie).zfill(6) + '.nc', 'w', format="NETCDF4") as file:
      nx = file.createDimension('nx', meps_avg.shape[1])
      ny = file.createDimension('ny', meps_avg.shape[0])
      nz = file.createDimension('nz_' + forecast_variable, 1)
      lat = file.createVariable('lat',np.float64,('ny','nx'))
      lon = file.createVariable('lon',np.float64,('ny','nx'))
      oro = file.createVariable('oro',np.float64,('ny','nx'))
      values = file.createVariable(forecast_variable,np.float64,('nz_' + forecast_variable,'ny','nx'))
      lat[:,:] = meps_lat[:,:]
      lon[:,:] = meps_lon[:,:]
      oro[:,:] = meps_oro[:,:]
      values[0,:,:] = meps_values[ie,:,:]
