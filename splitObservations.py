#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
import argparse
import numpy as np
import netCDF4
from datetime import datetime

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("full", help="Full obs.")
parser.add_argument("control", help="Control obs.")
parser.add_argument("assim", help="Assimilated obs.")
parser.add_argument("sigmao", help="Sigma obs")
args = parser.parse_args()

# Percentage of control observations
pC = 10.0

# Read all valid observations
with netCDF4.Dataset(args.full, 'r') as file:
  timeAll = file.variables['time'][:,:]
  locAll = file.variables['loc'][:,:]
  colsAll = file.variables['cols'][:,:]
  ncolAll = colsAll.shape[1]
  icolObs = -1
  icolErr = -1
  for icol in range(0, ncolAll):
    colName = getattr(file.variables['cols'], 'column_' + str(icol))
    if colName == "ObsVal":
      icolObs = icol
    if colName == "ObsErr":
      icolErr = icol
  if icolObs == -1:
    print("Cannot find ObsVal column")
    exit()
  if icolErr == -1:
    print("Cannot find ObsHofX column")
    exit()
obsSize = timeAll.shape[0]

# Number of control observations
nC = int(obsSize*pC/100.0)

# Random shuffle
arr = np.arange(obsSize)
np.random.shuffle(arr)

# Write control observations
with netCDF4.Dataset(args.control, 'w', format="NETCDF4") as file:
  nobs = file.createDimension('nobs', nC)
  ntime = file.createDimension('ntime', 6)
  nloc = file.createDimension('nloc', 3)
  ncol = file.createDimension('ncol', 1)
  time = file.createVariable('time',np.int32,('nobs', 'ntime'))
  loc = file.createVariable('loc',np.float64,('nobs', 'nloc'))
  cols = file.createVariable('cols',np.float64,('nobs', 'ncol'))
  cols.column_0 = "ObsVal"
  for jo in range(0, nC):
    io = arr[jo]
    time[jo,:] = timeAll[io,:]
    loc[jo,:] = locAll[io,:]
    cols[jo,0] = colsAll[io,icolObs]

# Write assimilated observations
with netCDF4.Dataset(args.assim, 'w', format="NETCDF4") as file:
  nobs = file.createDimension('nobs', obsSize-nC)
  ntime = file.createDimension('ntime', 6)
  nloc = file.createDimension('nloc', 3)
  ncol = file.createDimension('ncol', 2)
  time = file.createVariable('time',np.int32,('nobs', 'ntime'))
  loc = file.createVariable('loc',np.float64,('nobs', 'nloc'))
  cols = file.createVariable('cols',np.float64,('nobs', 'ncol'))
  cols.column_0 = "ObsVal"
  cols.column_1 = "ObsErr"
  for jo in range(nC, obsSize):
    io = arr[jo]
    time[jo-nC,:] = timeAll[io,:]
    loc[jo-nC,:] = locAll[io,:]
    cols[jo-nC,0] = colsAll[io,icolObs]
    cols[jo-nC,1] = float(args.sigmao)
