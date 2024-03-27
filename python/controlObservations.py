#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
  Prepare observations (screening and thinning):
   1. Read background
   2. Read observations
   3. Interpolate background in observation space
   4. Screen observations
   5. Thin observations
   6. Write observations
"""

#%%
import argparse
import gridpp
import numpy as np
import netCDF4
import titanlib
from datetime import datetime

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("background", help="Background file")
parser.add_argument("inputObs", help="Input observation file")
parser.add_argument("outputObs", help="Output observation file")
parser.add_argument("sigmab", help="Background error standard-deviation")
parser.add_argument("sigmao", help="Observation error standard-deviation")
parser.add_argument("scrFac", help="Screening factor")
args = parser.parse_args()

"""
   Parameters
"""

# Variable
forecast_variable = "air_temperature_2m"

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

# Read observations
with netCDF4.Dataset(args.inputObs, 'r') as file:
  timeAll = file.variables['time'][:,:]
  locAll = file.variables['loc'][:,:]
  obsAll = file.variables['cols'][:,0]
nobsAll = len(obsAll)

# Define observations points
points = gridpp.Points(locAll[:,1], locAll[:,0], locAll[:,2])

"""
   3. Interpolate background in observation space
"""

# Interpolate background in observation space
bkgAll = gridpp.bilinear(bgrid, points, bkg)

"""
   4. Screen observations
"""

# Get error standard-deviations
threshold = float(args.scrFac)*np.sqrt(float(args.sigmab)**2+float(args.sigmao)**2)

# Screen observations
obsScr = []
for jo in range(0, nobsAll):
  dep = obsAll[jo]-bkgAll[jo]
  if dep < threshold:
    obsScr.append(jo)
nobsScr = len(obsScr)

# Print number of screened observations
print("Screening: " + str(nobsScr) + " obs. kept out of " + str(nobsAll) + " (" + str(100.0*nobsScr/nobsAll) + " %)")

"""
   5. Thin observations
"""

# Define titanlib grid
net_points = titanlib.Points(locAll[obsScr,1], locAll[obsScr,0], locAll[obsScr,2])

# SCT settings
inner_radius = 50000
outer_radius = 100000
num_min = 5
num_max = 100
num_iterations = 1
num_min_prof = 20
dzmin = 100
dhmin = 10000
dz = 200
t2pos = np.full([nobsScr], 4)
t2neg = np.full([nobsScr], 4)
eps2 = np.full([nobsScr], 0.5)

# Run the SCT, method to exclude bad observations
flags, sct, rep = titanlib.sct(net_points, obsAll[obsScr],
        num_min, num_max, inner_radius, outer_radius, num_iterations,
        num_min_prof, dzmin, dhmin , dz, t2pos, t2neg, eps2)

# Get thinned observations
obsThn = []
for jo in range(0, nobsScr):
  if flags[jo] == 0:
    obsThn.append(obsScr[jo])
nobsThn = len(obsThn)

# Print number of thinned observations
print("Thinning:  " + str(nobsThn) + " obs. kept out of " + str(nobsScr) + " (" + str(100.0*nobsThn/nobsScr) + " %)")

"""
   6. Write observations
"""

# Write observations
with netCDF4.Dataset(args.outputObs, 'w', format="NETCDF4") as file:
  nobs = file.createDimension('nobs', nobsThn)
  ntime = file.createDimension('ntime', 6)
  nloc = file.createDimension('nloc', 3)
  ncolThn = file.createDimension('ncol', 3)
  timeThn = file.createVariable('time',np.int32,('nobs', 'ntime'))
  locThn = file.createVariable('loc',np.float64,('nobs', 'nloc'))
  colsThn = file.createVariable('cols',np.float64,('nobs', 'ncol'))
  colsThn.column_0 = "ObsVal"
  colsThn.column_1 = "ObsErr"
  colsThn.column_2 = "ObsBkg"
  for jo in range(0, nobsThn):
    io = obsThn[jo]
    timeThn[jo,:] = timeAll[io,:]
    locThn[jo,:] = locAll[io,:]
    colsThn[jo,0] = obsAll[io]
    colsThn[jo,1] = float(args.sigmao)
    colsThn[jo,2] = bkgAll[io]
