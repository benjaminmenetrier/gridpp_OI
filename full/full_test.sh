#!/bin/bash

# Number of processors
np=4

# Number of tests
nt=1

# Clean logs
rm -f full_log full_tests_results full_tests_diags

# Prepare background
#cp ../background_hr.nc full_background.nc

# Prepare observations
#cp ../observations.nc full_observations.nc

# Prepare background error
# - Standard-deviation
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_convertstate.x full_createStdDev.json >> full_log
# - BUMP
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_error_covariance_toolbox.x full_prepareBackgroundErrorBUMP.json >> full_log
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_error_covariance_toolbox.x full_diracBackgroundErrorBUMP.json >> full_log
# - FastLAM
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_error_covariance_toolbox.x full_prepareBackgroundErrorFastLAM.json >> full_log
#mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_error_covariance_toolbox.x full_diracBackgroundErrorFastLAM.json >> full_log

for it in $(seq 1 ${nt}); do
  # Split observation
#  python ../splitObservations.py full_observations.nc full_observations_control.nc full_observations_assim.nc 1.5 >> full_log

  # Background departure for control observations
#  ln -sf full_background.nc full_hofx_input_field.nc
#  ln -sf full_observations_control.nc full_hofx_input_obs.nc
#  mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_hofx.x full_hofx.json >> full_log
#  mv full_hofx_output_obs.nc full_observations_control_background.nc
#  python ../compareObservations.py background full_observations_control_background.nc >> full_tests_results

  # Run OI
#  python ../runOI.py full_background full_observations_assim full_analysis_OI >> full_log
#  ncdiff -O full_analysis_OI.nc full_background.nc full_increment_OI.nc >> full_log
#  ln -sf full_analysis_OI.nc full_hofx_input_field.nc
#  ln -sf full_observations_control.nc full_hofx_input_obs.nc
#  mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_hofx.x full_hofx.json >> full_log
#  mv full_hofx_output_obs.nc full_observations_control_OI.nc
#  python ../compareObservations.py OI full_observations_control_OI.nc >> full_tests_results

  # Screening
#  ln -sf full_background.nc full_hofx_input_field.nc
#  ln -sf full_observations_assim.nc full_hofx_input_obs.nc
#  mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_hofx.x full_hofx.json >> full_log
#  mv full_hofx_output_obs.nc full_observations_screening.nc
#  python ../screenObservations.py full_observations_screening.nc full_observations_screened.nc 2.7
  cp -f full_observations_assim.nc full_observations_screened.nc

  # Run 3DVar
  mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_variational.x full_3dvarFastLAM.json >> full_log
  ncdiff -O full_analysis_3dvar.nc full_background.nc full_increment_3dvar.nc >> full_log
  ln -sf full_analysis_3dvar.nc full_hofx_input_field.nc
  ln -sf full_observations_control.nc full_hofx_input_obs.nc
  mpirun -n ${np} ${HOME}/build/oops-bundle/bin/quench_hofx.x full_hofx.json >> full_log
  mv full_hofx_output_obs.nc full_observations_control_3dvar.nc
  python ../compareObservations.py 3dvar full_observations_control_3dvar.nc >> full_tests_results

  # Run Desroziers' diagnostics
#  python ../diagnoseDesroziers.py full_observations_assim_3dvar.nc >> full_tests_diags
done

# Compare OI/3dvar
python full_plot.py
