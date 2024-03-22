#!/bin/bash

# Number of tests
nT=1

# Cleaning
#rm -f log
#rm -f test_results

# Get data
#python getData.py >> log

# Background departure for all observations
#mpirun -n 4 ${HOME}/build/oops-bundle/bin/quench_hofx.x hofx_background.json >> log
#python compareObservations.py background observations_background.nc

# Prepare background error
#mpirun -n 4 ${HOME}/build/oops-bundle/bin/quench_convertstate.x createStdDev.json >> log
#mpirun -n 4 ${HOME}/build/oops-bundle/bin/quench_error_covariance_toolbox.x prepareBackgroundError.json
#mpirun -n 4 ${HOME}/build/oops-bundle/bin/quench_error_covariance_toolbox.x diracBackgroundError.json

for iT in $(seq 1 ${nT}); do
  # Split observation
  echo "Split observation - "${iT}
#  python splitObservations.py >> log

  # Background departure for control observations
  echo "Background departure for control observations - "${iT}
#  mpirun -n 4 ${HOME}/build/oops-bundle/bin/quench_hofx.x hofx_test_background.json >> log
#  python compareObservations.py background observations_control_background.nc >> test_results

  # Run OI
  echo "Run OI - "${iT}
#  python runOI.py background_hr observations_assim analysis_OI >> log
#  ncdiff -O analysis_OI.nc background_hr.nc increment_OI.nc >> log

  # OI analysis departure for control observations
  echo "OI analysis departure for control observations - "${iT}
#  mpirun -n 4 ${HOME}/build/oops-bundle/bin/quench_hofx.x hofx_test_OI.json >> log
#  python compareObservations.py OI observations_control_OI.nc >> test_results

  # Run 3DVar
  echo "Run 3DVar - "${iT}
  mpirun -n 4 ${HOME}/build/oops-bundle_debug/bin/quench_variational.x 3dvar_test.json
  ncdiff -O analysis_3dvar.nc background_hr.nc increment_3dvar.nc >> log

  # 3DVar analysis departure for control observations
  echo "3DVar analysis departure for control observations - "${iT}
#  mpirun -n 4 ${HOME}/build/oops-bundle/bin/quench_hofx.x hofx_test_3dvar.json >> log
#  python compareObservations.py 3DVar observations_control_3dvar.nc >> test_results
done
