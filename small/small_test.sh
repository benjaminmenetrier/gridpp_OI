#!/bin/bash

# Number of processors
np=4

# Number of tests
nt=1

# Clean logs
rm -f small_log small_tests_results small_tests_diags

# Prepare background
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_convertstate.x small_prepareBackground.json >> small_log

# Prepare observations
#python small_prepareObservations.py >> small_log

# Hollingworgh-Lönnberg diagnostics
#ln -sf small_background.nc small_hofx_input_field.nc
#ln -sf small_observations.nc small_hofx_input_obs.nc
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_hofx.x small_hofx.json >> small_log
#mv small_hofx_output_obs.nc small_observations_background.nc
#python ../diagnoseHL.py small_observations_background.nc

# Prepare background error
# - Standard-deviation
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_convertstate.x small_createStdDev.json >> small_log
# - BUMP
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_error_covariance_toolbox.x small_prepareBackgroundErrorBUMP.json >> small_log
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_error_covariance_toolbox.x small_diracBackgroundErrorBUMP.json >> small_log
# - FastLAM
mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_error_covariance_toolbox.x small_prepareBackgroundErrorFastLAM.json >> small_log
mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_error_covariance_toolbox.x small_diracBackgroundErrorFastLAM.json >> small_log

for it in $(seq 1 ${nt}); do
  # Split observation
  python ../splitObservations.py small_observations.nc small_observations_control.nc small_observations_assim.nc 1.5 >> small_log

  # Background departure for control observations
  ln -sf small_background.nc small_hofx_input_field.nc
  ln -sf small_observations_control.nc small_hofx_input_obs.nc
  mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_hofx.x small_hofx.json >> small_log
  mv small_hofx_output_obs.nc small_observations_control_background.nc
  python ../compareObservations.py background small_observations_control_background.nc >> small_tests_results

  # Run OI
  python ../runOI.py small_background small_observations_assim small_analysis_OI >> small_log
  ncdiff -O small_analysis_OI.nc small_background.nc small_increment_OI.nc >> small_log
  ln -sf small_analysis_OI.nc small_hofx_input_field.nc
  ln -sf small_observations_control.nc small_hofx_input_obs.nc
  mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_hofx.x small_hofx.json >> small_log
  mv small_hofx_output_obs.nc small_observations_control_OI.nc
  python ../compareObservations.py OI small_observations_control_OI.nc >> small_tests_results

  # Run 3DVar
  mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_variational.x small_3dvarFastLAM.json >> small_log
  ncdiff -O small_analysis_3dvar.nc small_background.nc small_increment_3dvar.nc >> small_log
  ln -sf small_analysis_3dvar.nc small_hofx_input_field.nc
  ln -sf small_observations_control.nc small_hofx_input_obs.nc
  mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_hofx.x small_hofx.json >> small_log
  mv small_hofx_output_obs.nc small_observations_control_3dvar.nc
  python ../compareObservations.py 3dvar small_observations_control_3dvar.nc >> small_tests_results

  # Run Desroziers' diagnostics
#  python ../diagnoseDesroziers.py small_observations_assim_3dvar.nc >> small_tests_diags
done

# Compare OI/3dvar
python small_plot.py
