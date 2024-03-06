#!/bin/bash

# Number of tests
nT=1

# Get data
#python getData.py

# Background departure for all observations
#mpirun -n 4 ~/build/oops-bundle/bin/quench_hofx.x hofx_background.json
#python compareObservations.py background observations_background.nc

# Split operations
#python splitObservations.py ${nT}

# Run over TEST
for iT in $(seq 1 ${nT}); do
  # Background departure for control observations
#  sed -e s/_TEST_/${iT}/g -e s/_INPUT_/background_hr/g -e s/_OUTPUT_/observations_control_background_${iT}/g hofx_template.json > hofx_test.json
#  mpirun -n 4 ~/build/oops-bundle/bin/quench_hofx.x hofx_test.json
#  python compareObservations.py background observations_control_background_${iT}.nc

  # Run 3DVar
#  sed -e s/_TEST_/${iT}/g 3dvar_template.json > 3dvar_test.json
#  mpirun -n 4 ~/build/oops-bundle/bin/quench_variational.x 3dvar_test.json
#  ncdiff -O analysis_3dvar_${iT}.nc background_hr.nc increment_3dvar_${iT}.nc

  # 3DVar analysis departure for control observations
#  sed -e s/_TEST_/${iT}/g -e s/_INPUT_/analysis_3dvar_${iT}/g -e s/_OUTPUT_/observations_control_3dvar_${iT}/g hofx_template.json > hofx_test.json
#  mpirun -n 4 ~/build/oops-bundle/bin/quench_hofx.x hofx_test.json
#  python compareObservations.py 3DVar observations_control_3dvar_${iT}.nc

  # Run OI
  python runOI.py observations_assim_${iT} analysis_OI_${iT}
  ncdiff -O analysis_OI_${iT}.nc background_hr.nc increment_OI_${iT}.nc

  # OI analysis departure for control observations
  sed -e s/_TEST_/${iT}/g -e s/_INPUT_/analysis_OI_${iT}/g -e s/_OUTPUT_/observations_control_OI_${iT}/g hofx_template.json > hofx_test.json
  mpirun -n 4 ~/build/oops-bundle/bin/quench_hofx.x hofx_test.json
  python compareObservations.py OI observations_control_OI_${iT}.nc
done
