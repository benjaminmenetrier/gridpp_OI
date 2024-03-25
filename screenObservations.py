#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

#%%
import argparse
import numpy as np
import netCDF4
from datetime import datetime

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input file")
parser.add_argument("output", help="Output file")
parser.add_argument("sigmab", help="Sigma b")
args = parser.parse_args()

# Read observations
with netCDF4.Dataset(args.input, 'r') as file:
  time = file.variables['time'][:,:]
  loc = file.variables['loc'][:,:]
  cols = file.variables['cols'][:,:]
  ncol = cols.shape[1]
  icolObs = -1
  icolErr = -1
  icolBkg = -1
  for icol in range(0, ncol):
    colName = getattr(file.variables['cols'], 'column_' + str(icol))
    if colName == "ObsVal":
      icolObs = icol
    if colName == "ObsErr":
      icolErr = icol
    if colName == "ObsHofX":
      icolBkg = icol
  if icolObs == -1:
    print("Cannot find ObsVal column")
    exit()
  if icolErr == -1:
    print("Cannot find ObsErr column")
    exit()
  if icolBkg == -1:
    print("Cannot find ObsHofX column")
    exit()

# Screen observations
sigmab = float(args.sigmab)
validObs = []
for jobs in range(0, cols.shape[0]):
  normDep = cols[jobs,icolObs]-cols[jobs,icolBkg]
  sigmao = cols[jobs,icolErr]
  if normDep < 1.0*(sigmab+sigmao):
    validObs.append(jobs)

# Print number of screened observations
print("Screening: " + str(len(validObs)) + " obs. kept out of " + str(cols.shape[0]))

# Write screened observations
with netCDF4.Dataset(args.output, 'w', format="NETCDF4") as file:
  nobs = file.createDimension('nobs', len(validObs))
  ntime = file.createDimension('ntime', 6)
  nloc = file.createDimension('nloc', 3)
  ncolScr = file.createDimension('ncol', 2)
  timeScr = file.createVariable('time',np.int32,('nobs', 'ntime'))
  locScr = file.createVariable('loc',np.float64,('nobs', 'nloc'))
  colsScr = file.createVariable('cols',np.float64,('nobs', 'ncol'))
  colsScr.column_0 = "ObsVal"
  colsScr.column_1 = "ObsErr"
  for jo in range(0, len(validObs)):
    io = validObs[jo]
    timeScr[jo,:] = time[io,:]
    locScr[jo,:] = loc[io,:]
    colsScr[jo,0] = cols[io,icolObs]
    colsScr[jo,1] = cols[io,icolErr]
