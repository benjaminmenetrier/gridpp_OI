#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
  Compute errors with control observations
   1. Read observations and background in observation space 
   2. Read background grid, OI and 3dvar analyses
   3. Interpolate OI and 3dvar analyses in observation space
   4. Compute and write errors
"""

import os
import argparse
import gridpp
import numpy as np
import netCDF4

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("obs", help="Observations")
parser.add_argument("background", help="Background")
parser.add_argument("analysisOI", help="OI analysis")
parser.add_argument("analysisVar", help="3dvar analysis")
parser.add_argument("errors", help="Errors file")
args = parser.parse_args()

"""
   Parameters
"""

# Variable
forecast_variable = "air_temperature_2m"

"""
   1. Read observations and background in observation space 
"""

# Read observations and background in observation space
with netCDF4.Dataset(args.obs, 'r') as file:
  loc = file.variables['loc'][:,:]
  cols = file.variables['cols'][:,:]
  ncol = cols.shape[1]
  icolObs = -1
  icolBkg = -1
  for icol in range(0, ncol):
    colName = getattr(file.variables['cols'], 'column_' + str(icol))
    if colName == "ObsVal":
      icolObs = icol
    if colName == "ObsBkg":
      icolBkg = icol
  if icolObs == -1:
    print("Cannot find ObsVal column")
    exit()
  if icolBkg == -1:
    print("Cannot find ObsHofX column")
    exit()
  obs = cols[:,icolObs]
  bkg = cols[:,icolBkg]
  lon = loc[:,0]
  lat = loc[:,1]
nobs = cols.shape[0]

"""
   2. Read background grid, OI and 3dvar analyses
"""

# Read background grid
with netCDF4.Dataset(args.background, 'r') as file:
  blat = file.variables['lat'][:,:]
  blon = file.variables['lon'][:,:]
  boro = file.variables['oro'][0,:,:]

# Read OI analysis
with netCDF4.Dataset(args.analysisOI, 'r') as file:
  ana_oi = file.variables[forecast_variable][0,:,:]

# Read 3dvar analysis
with netCDF4.Dataset(args.analysisVar, 'r') as file:
  ana_var = file.variables[forecast_variable][0,:,:]

"""
   3. Interpolate OI and 3dvar analyses in observation space
"""

# Define observations points
points = gridpp.Points(loc[:,1], loc[:,0], loc[:,2])

# Define grid
bgrid = gridpp.Grid(blat, blon, boro)

# Interpolate OI and 3dvar analyses in observation space
oi = gridpp.bilinear(bgrid, points, ana_oi)
var = gridpp.bilinear(bgrid, points, ana_var)

"""
   4. Compute and write errors
"""

# Compute error
bkgErr = np.abs(bkg-obs)
oiErr = np.abs(oi-obs)
varErr = np.abs(var-obs)
diffErr = varErr-oiErr
diffErrNeg = np.count_nonzero(diffErr<0)/len(diffErr)

# Print statistics
print("Bkg:  " + str(np.mean(bkgErr)) + " +/- " + str(np.std(bkgErr)))
print("OI:   " + str(np.mean(oiErr)) + " +/- " + str(np.std(oiErr)))
print("Var:  " + str(np.mean(varErr)) + " +/- " + str(np.std(varErr)))
print("Diff: " + str(np.mean(diffErr)) + " +/- " + str(np.std(diffErr)) + " (" + str(100.0*diffErrNeg) + " %)")

# Write errors
with netCDF4.Dataset(args.errors, 'w', format="NETCDF4") as file:
  n = file.createDimension('n', None)
  lonNC = file.createVariable('lon',np.float64,('n'))
  latNC = file.createVariable('lat',np.float64,('n'))
  bkgErrNC = file.createVariable('bkg',np.float64,('n'))
  oiErrNC = file.createVariable('oi',np.float64,('n'))
  varErrNC = file.createVariable('var',np.float64,('n'))
  lonNC[:] = lon[:]
  latNC[:] = lat[:]
  bkgErrNC[:] = bkgErr[:]
  oiErrNC[:] = oiErr[:]
  varErrNC[:] = varErr[:]
