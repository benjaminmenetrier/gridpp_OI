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
parser.add_argument("nT", help="Number of tests")

# Parse and print arguments
args = parser.parse_args()
print("Parameters:")
for arg in vars(args):
    if not arg is None:
        print(" - " + arg + ": " + str(getattr(args, arg)))

# Percentage of control observations
pC = 10.0

# Read all valid observations
with netCDF4.Dataset('observations.nc', 'r') as file:
    timeAll = file.variables['time'][:,:]
    locAll = file.variables['loc'][:,:]
    colsAll = file.variables['cols'][:,:]
obsSize = timeAll.shape[0]

# Number of control observations
nC = int(obsSize*pC/100.0)

for iT in range(0, int(args.nT)):
    # Random shuffle
    arr = np.arange(obsSize)
    np.random.shuffle(arr)

    # Write control observations
    with netCDF4.Dataset('observations_control_' + str(iT+1) + '.nc', 'w', format="NETCDF4") as file:
        nobs = file.createDimension('nobs', nC)
        ntime = file.createDimension('ntime', 6)
        nloc = file.createDimension('nloc', 3)
        ncol = file.createDimension('ncol', 1)
        time = file.createVariable('time',np.int32,('nobs', 'ntime'))
        loc = file.createVariable('loc',np.float64,('nobs', 'nloc'))
        cols = file.createVariable('cols',np.float64,('nobs', 'ncol'))
        cols.column_0 = "ObsVal"
        for jo in range(0, nC):
            jos = arr[jo]
            time[jo,:] = timeAll[jos,:]
            loc[jo,:] = locAll[jos,:]
            cols[jo,:] = colsAll[jos,0]

    # Write assimilated observations
    with netCDF4.Dataset('observations_assim_' + str(iT+1) + '.nc', 'w', format="NETCDF4") as file:
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
            jos = arr[jo]
            time[jo-nC,:] = timeAll[jos,:]
            loc[jo-nC,:] = locAll[jos,:]
            cols[jo-nC,:] = colsAll[jos,:]
