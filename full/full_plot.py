#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import numpy as np

# Prepare arrays
oi = []
var = []

# Read results
with open('full_tests_results', 'r') as f:
  lines = [line.rstrip() for line in f]
  for line in lines:
    test = line.split()[0]
    value = float(line.split()[1])
    if test == "background:":
      bkg = value
    elif test == "OI:":
      oi.append(value/bkg)
    elif test == "3dvar:":
      var.append(value/bkg)

# Compute statistics
oi_np = np.array(oi)
oi_mean = np.average(oi_np)
oi_std = np.std(oi_np)
var_np = np.array(var)
var_mean = np.average(var_np)
var_std = np.std(var_np)

# Print results
print("OI:    " + str(oi_mean) + " +/- " + str(oi_std))
print("3dvar: " + str(var_mean) + " +/- " + str(var_std))
