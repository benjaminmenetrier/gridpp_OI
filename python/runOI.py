#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
  Run gridpp-based OI:
   1. Read background
   2. Read observations
   3. Interpolate background in observation space
   4. Prepare structure functions
   5. Run OI
   6. Write analysis
"""

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
parser.add_argument("analysisInObsSpace", help="Analysis in obs space file")
args = parser.parse_args()

"""
   Parameters
"""

# Variable
forecast_variable = "air_temperature_2m"

# Variance ratio (simgao**2/sigmab**2)
variance_ratio = 0.1

# Correlation in horizontal and vertical directions
h = 25000
v = 500

"""
   1. Read background
"""

# Read background
with netCDF4.Dataset(args.background, 'r') as file:
  blat = file.variables['lat'][:,:]
  blon = file.variables['lon'][:,:]
  boro = file.variables['oro'][0,:,:]
  bkg = file.variables[forecast_variable][0,:,:]

# Define background grid
bgrid = gridpp.Grid(blat, blon, boro)

"""
   2. Read observations
"""

# Load observations
with netCDF4.Dataset(args.observations, 'r') as file:
  timeAll = file.variables['time'][:,:]
  locAll = file.variables['loc'][:,:]
  colsAll = file.variables['cols'][:,:]
  ncolAll = colsAll.shape[1]
  colName = {}
  colIndex = {}
  for icol in range(0, ncolAll):
    colAttr = 'column_' + str(icol)
    colName[colAttr] = getattr(file.variables['cols'], colAttr)
    colIndex[colName[colAttr]] = icol
nobsAll = colsAll.shape[0]

# Define observations points
points = gridpp.Points(locAll[:,1], locAll[:,0], locAll[:,2])

"""
   3. Interpolate background in observation space
"""

# Interpolate background in observation space
pbkg = gridpp.bilinear(bgrid, points, bkg)

"""
   4. Prepare structure functions
"""

# Variance ratios
variance_ratios = np.full(points.size(), variance_ratio)

# Prepare structuer functions
structure = gridpp.BarnesStructure(h, v)

"""
   5. Run OI
"""

# Optimal interpolation
max_points = 10
obs = colsAll[:,colIndex["ObsVal"]]
ana = gridpp.optimal_interpolation(bgrid, bkg, points, obs, variance_ratios, pbkg, structure, max_points)

"""
   6. Write analysis
"""

# Write analysis
with netCDF4.Dataset(args.analysis, 'w', format="NETCDF4") as file:
  nx = file.createDimension('nx', bkg.shape[1])
  ny = file.createDimension('ny', bkg.shape[0])
  nz = file.createDimension('nz_' + forecast_variable, 1)
  lat = file.createVariable('lat',np.float64,('ny','nx'))
  lon = file.createVariable('lon',np.float64,('ny','nx'))
  values = file.createVariable(forecast_variable,np.float64,('nz_' + forecast_variable,'ny','nx'))
  lat[:,:] = blat[:,:]
  lon[:,:] = blon[:,:]
  values[0,:,:] = ana[:,:]
