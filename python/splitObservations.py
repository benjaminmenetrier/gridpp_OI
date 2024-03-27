#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
  Split observations between control and assimilated
   1. Read observations
   2. Draw control and assimilated subsets
   3. Write control observations
   4. Write assimilated observations
"""

import argparse
import numpy as np
import netCDF4
from datetime import datetime

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("all", help="All observations")
parser.add_argument("control", help="Control observations")
parser.add_argument("assimilated", help="Assimilated observations")
args = parser.parse_args()

"""
   Parameters
"""

# Percentage of control observations
pC = 10.0

"""
   1. Read observations
"""

# Read observations
with netCDF4.Dataset(args.all, 'r') as file:
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
   2. Draw control and assimilated subsets
"""

# Number of control observations
nobsCtrl = int(nobsAll*pC/100.0)
nobsAssim = nobsAll-nobsCtrl

# Random shuffle
arr = np.arange(nobsAll)
np.random.shuffle(arr)

# 
timeCtrl = timeAll[arr[0:nobsCtrl],:]
locCtrl = locAll[arr[0:nobsCtrl],:]
colsCtrl = colsAll[arr[0:nobsCtrl],:]
timeAssim = timeAll[arr[nobsCtrl:nobsAll],:]
locAssim = locAll[arr[nobsCtrl:nobsAll],:]
colsAssim = colsAll[arr[nobsCtrl:nobsAll],:]

# Print repartition of observations
print("Control:      " + str(nobsCtrl) + " obs. out of " + str(nobsAll) + " (" + str(100.0*nobsCtrl/nobsAll) + " %)")
print("Assimilation: " + str(nobsAssim) + " obs. out of " + str(nobsAll) + " (" + str(100.0*nobsAssim/nobsAll) + " %)")

"""
   3. Write control observations
"""

# Write control observations
with netCDF4.Dataset(args.control, 'w', format="NETCDF4") as file:
  nobs = file.createDimension('nobs', nobsCtrl)
  ntime = file.createDimension('ntime', 6)
  nloc = file.createDimension('nloc', 3)
  ncol = file.createDimension('ncol', ncolAll)
  time = file.createVariable('time',np.int32,('nobs', 'ntime'))
  loc = file.createVariable('loc',np.float64,('nobs', 'nloc'))
  cols = file.createVariable('cols',np.float64,('nobs', 'ncol'))
  for icol in range(0, ncolAll):
    colAttr = "column_" + str(icol)
    setattr(file.variables['cols'], colAttr, colName[colAttr])
  time[:,:] = timeCtrl[:,:]
  loc[:,:] = locCtrl[:,:]
  cols[:,:] = colsCtrl[:,:]

"""
   4. Write assimilated observations
"""

# Write assimilated observations
with netCDF4.Dataset(args.assimilated, 'w', format="NETCDF4") as file:
  nobs = file.createDimension('nobs', nobsAll-nobsCtrl)
  ntime = file.createDimension('ntime', 6)
  nloc = file.createDimension('nloc', 3)
  ncol = file.createDimension('ncol', ncolAll)
  time = file.createVariable('time',np.int32,('nobs', 'ntime'))
  loc = file.createVariable('loc',np.float64,('nobs', 'nloc'))
  cols = file.createVariable('cols',np.float64,('nobs', 'ncol'))
  for icol in range(0, ncolAll):
    colAttr = "column_" + str(icol)
    setattr(file.variables['cols'], colAttr, colName[colAttr])
  time[:,:] = timeAssim[:,:]
  loc[:,:] = locAssim[:,:]
  cols[:,:] = colsAssim[:,:]
