#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

#%%
import argparse
import gridpp
import numpy as np
import netCDF4
import datetime

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("background", help="Background file")
parser.add_argument("observations", help="Observations file")
parser.add_argument("analysis", help="Analysis file")
args = parser.parse_args()

# Variable
forecast_variable = "air_temperature_2m"

# Load background
with netCDF4.Dataset(args.background + '.nc', 'r') as file:
  blat = file.variables['lat'][:,:]
  blon = file.variables['lon'][:,:]
  if len(file.variables['oro'].shape) == 2:
    boro = file.variables['oro'][:,:]
  elif len(file.variables['oro'].shape) == 3:
    boro = file.variables['oro'][0,:,:]
  background = file.variables[forecast_variable][0,:,:]

# Define background grid
bgrid = gridpp.Grid(blat, blon, boro)

# Load observations
with netCDF4.Dataset(args.observations + '.nc', 'r') as file:
  loc = file.variables['loc'][:,:]
  obs_values = file.variables['cols'][:,0]

# Define observations points
points = gridpp.Points(loc[:,1], loc[:,0], loc[:,2])

# Get background in observation space
pbackground = gridpp.bilinear(bgrid, points, background)

if False:
  # Write background in observations space
  obsSize = pbackground.shape[0]
  with netCDF4.Dataset('observations_in_obs_space.nc', 'w', format="NETCDF4") as file:
    nobs = file.createDimension('nobs', obsSize)
    nloc = file.createDimension('nloc', 3)
    ncol = file.createDimension('ncol', 2)
    locBkg = file.createVariable('loc',np.float64,('nobs', 'nloc'))
    cols = file.createVariable('cols',np.float64,('nobs', 'ncol'))
    cols.column_0 = "ObsVal"
    cols.column_1 = "ObsHofX"
    for jo in range(0, obsSize):
      locBkg[jo,:] = loc[jo,:]
      cols[jo,0] = pbackground[jo]
      cols[jo,1] = obs_values[jo]

# We trust observations 10 times more than MEPS-forecast
variance_ratios = np.full(points.size(), 0.1)

# Correlation in xy and z-direction
h = 25000
v = 500
structure = gridpp.BarnesStructure(h, v)
max_points = 10

# Optimal interpolation
analysis = gridpp.optimal_interpolation(bgrid, background, points, obs_values,
           variance_ratios, pbackground, structure, max_points)

# Write analysis
with netCDF4.Dataset(args.analysis + '.nc', 'w', format="NETCDF4") as file:
  nx = file.createDimension('nx', background.shape[1])
  ny = file.createDimension('ny', background.shape[0])
  nz = file.createDimension('nz_' + forecast_variable, 1)
  lat = file.createVariable('lat',np.float64,('ny','nx'))
  lon = file.createVariable('lon',np.float64,('ny','nx'))
  oro = file.createVariable('oro',np.float64,('ny','nx'))
  values = file.createVariable(forecast_variable,np.float64,('nz_' + forecast_variable,'ny','nx'))
  lat[:,:] = blat[:,:]
  lon[:,:] = blon[:,:]
  oro[:,:] = boro[:,:]
  values[0,:,:] = analysis[:,:]
