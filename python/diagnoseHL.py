#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

#%%
from alive_progress import alive_bar
import argparse
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from scipy.spatial import KDTree
from sklearn.metrics import pairwise_distances

# Parser
parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input file")
args = parser.parse_args()

# Parameters
diagnostic_points = 1
diagnostic_radius = 50.0e3
max_distance = 40.0e3
distance_class_width = 2.0e3
max_height = 100.0
height_class_width = 100.0
minSampleSize = 1000
nAmp = 100
nRad = 100
rad = np.linspace(0.0, 3.0*max_distance, nRad)
req = 6371.0e3

# Classes setup
n_distance = int(max_distance/distance_class_width)
distance_classes = np.zeros((n_distance))
distance_classes[1:] = np.linspace(0, (n_distance-2)*distance_class_width, n_distance-1)+distance_class_width/2
n_height = int(max_height/height_class_width)
height_classes = np.zeros((n_height))
height_classes[1:] = np.linspace(0, (n_height-2)*height_class_width, n_height-1)+height_class_width/2

# Gaussian function for fit
def func(x, a, sigma): 
  return a*np.exp(-(x**2)/(2.0*sigma**2))

# Load background and observations
with netCDF4.Dataset(args.input, 'r') as file:
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
    print("Cannot find ObsBkg column")
    exit()
  obs_values = cols[:,icolObs]
  bkg_values = cols[:,icolBkg]

# Compute obs-bkg
omb = obs_values-bkg_values

# Number of observations
obsSize = len(obs_values)

# Convert coordinates to radians
for iobs in range(0, obsSize):
  loc[iobs,0] = np.radians(loc[iobs,0])
  loc[iobs,1] = np.radians(loc[iobs,1])

# Build vector of cartesian coordinates
cc = np.zeros((obsSize, 3))
for iobs in range(0, obsSize):
  cc[iobs,0] = (req+loc[iobs,2])*np.cos(loc[iobs,0])*np.cos(loc[iobs,1])
  cc[iobs,1] = (req+loc[iobs,2])*np.sin(loc[iobs,0])*np.cos(loc[iobs,1])
  cc[iobs,2] = (req+loc[iobs,2])*np.sin(loc[iobs,1])

# Build KDTree
T = KDTree(cc)

# Define a subsample of diagnostic points
arr = np.arange(obsSize)
np.random.shuffle(arr)

# Initialize loop
foundPoints = 0
iobs = 0

while foundPoints < diagnostic_points:
  # Find surrounding points
  lidx = T.query_ball_point(cc[arr[iobs],:],r=diagnostic_radius)

  if len(lidx) >= minSampleSize:
    # Initialize arrays
    cov = np.zeros((n_distance, n_height))
    mean1 = np.zeros((n_distance, n_height))
    mean2 = np.zeros((n_distance, n_height))
    card = np.zeros((n_distance, n_height))

    # Compute moments
    with alive_bar(len(lidx)) as bar:
      for jobs in range(0, len(lidx)):
        centerLatLon = [[loc[lidx[jobs],1], loc[lidx[jobs],0]]]
        cidx = T.query_ball_point(cc[lidx[jobs],:],r=max_distance)
        latLon = np.zeros((len(cidx), 2))
        for kobs in range(0, len(cidx)):
          latLon[kobs,:] = [loc[cidx[kobs],1], loc[cidx[kobs],0]]
        hdist = pairwise_distances(latLon, centerLatLon, metric='haversine')
        hdist = hdist*req
        for kobs in range(0, len(cidx)):
          if hdist[kobs,0] < 1.0e-12:
            i = 0
          else:
            i = int(hdist[kobs,0]/distance_class_width)+1
          vdist = np.abs(loc[cidx[kobs],2]-loc[lidx[jobs],2])
          if vdist < 1.0e-12:
            k = 0
          else:
            if n_height == 1:
              k = int(vdist/height_class_width)
            else:
              k = int(vdist/height_class_width)+1
          if i < n_distance and k < n_height:
            if card[i,k] > 0.0:
              cov[i,k] += card[i,k]/(card[i,k]+1.0)*(omb[lidx[jobs]]-mean1[i,k])*(omb[cidx[kobs]]-mean2[i,k])
            mean1[i,k] += 1.0/(card[i,k]+1.0)*(omb[lidx[jobs]]-mean1[i,k])
            mean2[i,k] += 1.0/(card[i,k]+1.0)*(omb[cidx[kobs]]-mean2[i,k])
            card[i,k] += 1.0
        # Update progression
        bar()

    # Normalize covariance
    for i in range(0, n_distance):
      for k in range(0, n_height):
        if card[i,k] > 1.0:
          cov[i,k] /= card[i,k]-1.0

    # Brute-force fit
    amp = np.linspace(0.0, cov[0,0], nAmp)
    mseMin = 1.0e16
    for iAmp in range(0, nAmp):
      for iRad in range(0, nRad):
        if rad[iRad] > 0.0:
          mse = 0.0
          for ic in range(1, n_distance):
            mse += (func(distance_classes[ic], amp[iAmp], rad[iRad])-cov[ic,0])**2
          if mse < mseMin:
            mseMin = mse
            ampOpt = amp[iAmp]
            radOpt = rad[iRad]

    # Compute parameters
    sigmao = np.sqrt(cov[0,0]-ampOpt)
    sigmab = np.sqrt(ampOpt)
    rh = radOpt*3.57
    print("sigma_o: " + str(sigmao))
    print("sigma_b: " + str(sigmab))
    print("r_h: " + str(rh))
 
    # Done! 
    foundPoints += 1

  # Update diagnostic point index
  iobs += 1

exit()
# Curve fit
cov_fit = func(distance_classes, ampOpt, radOpt)

# Plot raw data
fig,ax = plt.subplots(ncols=1, figsize=(6,4))
ax.set_xlabel("Distance")
ax.set_ylabel("HBHT + R")
ax.set_xlim([0,1.1*distance_classes[n_distance-1]])
ax.set_ylim([0,1.1*cov[0,0]])
ax.scatter(distance_classes, cov[:,0])
ax.plot(distance_classes, cov_fit, c='r')
plt.show()
#fig.savefig('diagnoseHL.png')
