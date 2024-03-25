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
parser.add_argument("prefix", help="Standard output prefix")
parser.add_argument("file", help="HofX file")
args = parser.parse_args()

# Read observations
with netCDF4.Dataset(args.file, 'r') as file:
  cols = file.variables['cols'][:,:]
  ncol = cols.shape[1]
  icolTest = -1
  icolCtrl = -1
  for icol in range(0, ncol):
    colName = getattr(file.variables['cols'], 'column_' + str(icol))
    if colName == "ObsHofX":
      icolTest = icol
    if colName == "ObsVal":
      icolCtrl = icol
  if icolTest == -1:
    print("Cannot find ObsHofX column")
    exit()
  if icolCtrl == -1:
    print("Cannot find ObsVal column")
    exit()

# Compute RMSE
rmse = np.sqrt(np.average(np.square(cols[:,icolTest]-cols[:,icolCtrl])))
print(args.prefix + ": " + str(rmse))
