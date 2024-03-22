#!/bin/bash

# Number of processors
np=4

# Prepare background
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_convertstate.x small_prepareBackground.json 

# Prepare observations
python small_prepareObservations.py 

# Run OI
# - Run gridpp
python runOI.py small_background small_observations small_analysis_OI
# - Compute increment
ncdiff -O small_analysis_OI.nc small_background.nc small_increment_OI.nc

# Prepare background error
# - Standard-deviation
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_convertstate.x small_createStdDev.json 
# - BUMP
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_error_covariance_toolbox.x small_prepareBackgroundErrorBUMP.json
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_error_covariance_toolbox.x small_prepareBackgroundErrorBUMP.json
# - FastLAM
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_error_covariance_toolbox.x small_diracBackgroundErrorFastLAM.json
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_error_covariance_toolbox.x small_diracBackgroundErrorFastLAM.json 

# Run hofx
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_hofx.x small_hofx.json

# Run 3DVar
# - BUMP
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_variational.x small_3dvarBUMP.json
# - FastLAM
mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_variational.x small_3dvarFastLAM.json
# - Compute increment
ncdiff -O small_analysis_3dvar.nc small_background.nc small_increment_3dvar.nc

# Compare
ncdiff -O small_increment_3dvar.nc small_increment_OI.nc small_increment_diff.nc
