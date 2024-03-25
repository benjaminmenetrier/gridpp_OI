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
diagnostic_points = 1000
diagnostic_radius = 20.0e3
max_distance = 20.0e3
distance_class_width = 2.0e3
max_height = 100.0
height_class_width = 100.0
req = 6371.0e3

# Classes setup
n_distance = int(max_distance/distance_class_width)
distance_classes = np.zeros((n_distance))
distance_classes[1:] = np.linspace(0, (n_distance-2)*distance_class_width, n_distance-1)+distance_class_width/2
n_height = int(max_height/height_class_width)
height_classes = np.zeros((n_height))
height_classes[1:] = np.linspace(0, (n_height-2)*height_class_width, n_height-1)+height_class_width/2

# Load background
with netCDF4.Dataset(args.input, 'r') as file:
  loc = file.variables['loc'][:,:]
  cols = file.variables['cols'][:,:]
  ncol = cols.shape[1]
  icolOmbg = -1
  icolOman = -1
  for icol in range(0, ncol):
    colName = getattr(file.variables['cols'], 'column_' + str(icol))
    if colName == "ombg":
      icolOmbg = icol
    if colName == "oman":
      icolOman = icol
  if icolOmbg == -1:
    print("Cannot find ombg column")
    exit()
  if icolOman == -1:
    print("Cannot find oman column")
    exit()
  ombg = cols[:,icolOmbg]
  oman = cols[:,icolOman]
anmbg = ombg-oman

# Number of observations
obsSize = len(ombg)

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
didx = arr[0:diagnostic_points]

# Define local trees
lidx = []
totalSize = 0
for iobs in range(0, diagnostic_points):
  # Find all observations in the diagnostic radius
  lidx.append(T.query_ball_point(cc[didx[iobs],:],r=diagnostic_radius))
  totalSize += len(lidx[iobs])

# Initialize arrays
var = np.zeros((diagnostic_points))
meanOman = np.zeros((diagnostic_points))
cov = np.zeros((n_distance, n_height, diagnostic_points))
meanAnmbg = np.zeros((n_distance, n_height, diagnostic_points))
meanOmbg = np.zeros((n_distance, n_height, diagnostic_points))
card = np.zeros((n_distance, n_height, diagnostic_points))

with alive_bar(totalSize) as bar:
  for iobs in range(0, diagnostic_points):
    # Compute moments
    for jobs in range(0, len(lidx[iobs])):
      centerLatLon = [[loc[lidx[iobs][jobs],1], loc[lidx[iobs][jobs],0]]]
      cidx = T.query_ball_point(cc[lidx[iobs][jobs],:],r=max_distance)
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
        vdist = np.abs(loc[cidx[kobs],2]-loc[lidx[iobs][jobs],2])
        if vdist < 1.0e-12:
          k = 0
        else:
          k = int(vdist/height_class_width)+1       
        if i == 0 and k == 0:
          if card[0,0,iobs] > 0.0:
            var[iobs] += card[0,0,iobs]/(card[0,0,iobs]+1.0)*(oman[cidx[kobs]]-meanOman[iobs])*(ombg[lidx[iobs][jobs]]-meanOmbg[0,0,iobs])
          meanOman[iobs] += 1.0/(card[0,0,iobs]+1.0)*(oman[cidx[kobs]]-meanOman[iobs])
        if i < n_distance and k < n_height:
          if card[i,k,iobs] > 0.0:
            cov[i,k,iobs] += card[i,k,iobs]/(card[i,k,iobs]+1.0)*(anmbg[cidx[kobs]]-meanAnmbg[i,k,iobs])*(ombg[lidx[iobs][jobs]]-meanOmbg[i,k,iobs])
          meanAnmbg[i,k,iobs] += 1.0/(card[i,k,iobs]+1.0)*(anmbg[cidx[kobs]]-meanAnmbg[i,k,iobs])
          meanOmbg[i,k,iobs] += 1.0/(card[i,k,iobs]+1.0)*(ombg[lidx[iobs][jobs]]-meanOmbg[i,k,iobs])
          card[i,k,iobs] += 1.0
      bar()

# Average covariances
var_avg = 0.0
cov_avg = np.zeros((n_distance, n_height))
card_avg = np.zeros((n_distance, n_height))
for iobs in range(0, diagnostic_points):
  var_avg += var[iobs]
  for i in range(0, n_distance):
    for k in range(0, n_height):
      cov_avg[i,k] += cov[i,k,iobs]
      card_avg[i,k] += card[i,k,iobs]

# Compute covariance
if card_avg[0,0] > 1.0:
  var_avg /= card_avg[0,0]-1.0
print("var_avg: " + str(var_avg))
for i in range(0, n_distance):
  for k in range(0, n_height):
    if card_avg[i,k] > 1.0:
      cov_avg[i,k] /= card_avg[i,k]-1.0
print("cov_avg: " + str(cov_avg[:,0]))

# Gaussian function for fit
def func(x, a, sigma): 
  return a*np.exp(-(x**2)/(2.0*sigma**2))

# Brute-force fit
nAmp = 100
amp = np.linspace(0.0, cov_avg[0], nAmp)
nRad = 100
rad = np.linspace(0.0, 3.0*max_distance, nRad)
mseMin = 1.0e16
for iAmp in range(0, nAmp):
  for iRad in range(0, nRad):
    if rad[iRad] > 0.0:
      mse = 0.0
      for ic in range(0, n_distance):
        mse += (func(distance_classes[ic], amp[iAmp], rad[iRad])-cov_avg[ic,0])**2
      if mse < mseMin:
        mseMin = mse
        ampOpt = amp[iAmp]
        radOpt = rad[iRad]

# Curve fit
cov_fit = func(distance_classes, ampOpt, radOpt)

# Print results
print("Sigma_o: " + str(np.sqrt(var_avg)))
print("Sigma_b: " + str(np.sqrt(ampOpt)))
print("L_b: " + str(radOpt*3.57))

# Plot raw data
fig,ax = plt.subplots(ncols=1, figsize=(6,4))
ax.set_xlabel("Distance")
ax.set_ylabel("HBHT")
ax.set_xlim([0,1.1*distance_classes[n_distance-1]])
ax.set_ylim([0,1.1*cov_avg[0,0]])
ax.scatter(distance_classes, cov_avg[:,0])
ax.plot(distance_classes, cov_fit, c='r')
plt.show()
#fig.savefig('diagnoseDesroziers.png')
