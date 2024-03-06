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

# Variable
forecast_variable = "air_temperature_2m"

# Load background
with netCDF4.Dataset('background_hr.nc', 'r') as file:
    blats = file.variables['lats'][:,:]
    blons = file.variables['lons'][:,:]
    belevs = file.variables['elevs'][:,:]
    background = file.variables[forecast_variable][0,:,:]

# Define background grid
bgrid = gridpp.Grid(blats, blons, belevs)

# Single observation location
loc = np.array([[2.0, 70.0, 0.0]])

# Define observations points
points = gridpp.Points(loc[:,1], loc[:,0], loc[:,2])

# Get background in observation space
pbackground = gridpp.bilinear(bgrid, points, background)

# Single observation increment
obs_values = np.array([pbackground[0]+2.0])

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
variance_ratios = np.full(points.size(), 1.0)

# Correlation in xy and z-direction
h = 25000
v = 500
structure = gridpp.BarnesStructure(h, v)
max_points = 10

# Optimal interpolation
analysis = gridpp.optimal_interpolation(bgrid, background, points, obs_values,
           variance_ratios, pbackground, structure, max_points)

# Write analysis
with netCDF4.Dataset('single_obs_OI.nc', 'w', format="NETCDF4") as file:
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
    values[0,:,:] = analysis[:,:]-background[:,:]
