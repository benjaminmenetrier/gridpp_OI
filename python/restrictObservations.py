#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 Restrict observations to small domain:
   1. Read background coordinates and get domain bounds
   2. Read observations
   3. Select observations inside the domain
   4. Write observations
"""

import argparse
import numpy as np
import netCDF4
from datetime import datetime

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("background", help="Background file")
parser.add_argument("inputObs", help="Input observations file")
parser.add_argument("outputObs", help="Output observations file")
args = parser.parse_args()

"""
   Parameters
"""

# Variable
forecast_variable = "air_temperature_2m"

"""
   1. Read background coordinates and get domain bounds
"""

# Read background
with netCDF4.Dataset(args.background, 'r') as file:
  blat = file.variables['lat'][:,:]
  blon = file.variables['lon'][:,:]
nx = blon.shape[1]
ny = blon.shape[0]

# Get domain bounds
lonMinVec = []
lonMaxVec = []
for ix in range(0, nx):
  lonMinVec.append(np.min(blon[:,ix]))
  lonMaxVec.append(np.max(blon[:,ix]))
lonMin = min(lonMaxVec)
lonMax = max(lonMinVec)
latMinVec = []
latMaxVec = []
for iy in range(0, ny):
  latMinVec.append(np.min(blat[iy,:]))
  latMaxVec.append(np.max(blat[iy,:]))
latMin = min(latMaxVec)
latMax = max(latMinVec)

"""
   2. Read observations
"""

# Read observations
with netCDF4.Dataset(args.inputObs, 'r') as file:
  timeAll = file.variables['time'][:,:]
  locAll = file.variables['loc'][:,:]
  colsAll = file.variables['cols'][:,:]
  ncolAll = colsAll.shape[1]
  colName = {}
  for icol in range(0, ncolAll):
    colAttr = 'column_' + str(icol)
    colName[colAttr] = getattr(file.variables['cols'], colAttr)
nobsAll = colsAll.shape[0]

"""
   3. Select observations inside the domain
"""

# Keep observations inside the domain
obsIn = []
for jo in range(0, locAll.shape[0]):
  if (lonMin <= locAll[jo,0]) and (locAll[jo,0] <= lonMax) and (latMin <= locAll[jo,1]) and (locAll[jo,1] <= latMax):
    obsIn.append(jo)
nobsIn = len(obsIn)

# Print number of inside observation
print("Inside:  " + str(nobsIn) + " obs. kept out of " + str(nobsAll) + " (" + str(100.0*nobsIn/nobsAll) + " %)")

"""
   4. Write observations
"""

# Write observations
with netCDF4.Dataset(args.outputObs, 'w', format="NETCDF4") as file:
  nobs = file.createDimension('nobs', nobsIn)
  ntime = file.createDimension('ntime', 6)
  nloc = file.createDimension('nloc', 3)
  ncol = file.createDimension('ncol', ncolAll)
  time = file.createVariable('time',np.int32,('nobs', 'ntime'))
  loc = file.createVariable('loc',np.float64,('nobs', 'nloc'))
  cols = file.createVariable('cols',np.float64,('nobs', 'ncol'))
  for icol in range(0, ncolAll):
    colAttr = "column_" + str(icol)
    setattr(file.variables['cols'], colAttr, colName[colAttr])
  for jo in range(0, nobsIn):
    io = obsIn[jo]
    time[jo,:] = timeAll[io,:]
    loc[jo,:] = locAll[io,:]
    cols[jo,:] = colsAll[io,:]
