#!/bin/bash

update_json () {
  sed -e s/_PREFIX_/${prefix}/g \
      -e s:_DATADIR_:${datadir}:g \
      -e s/_NX_/${nx}/g \
      -e s/_NY_/${ny}/g \
      -e s/_LONCENTRE_/${lonCentre}/g \
      -e s/_LATCENTRE_/${latCentre}/g \
      -e s/_LON0_/${lon0}/g \
      -e s/_LAT0_/${lat0}/g \
      -e s/_SIGMAB_/${sigmab}/g \
      -e s/_RH_/${rh}/g \
      -e s/_RV_/${rv}/g \
      template/$1 > ${prefix}/$1
}

# Number of processors
np=8

# Number of tests
nt=5

# Quench binaries directory
bindir=${HOME}/build/oops-bundle/bin

# Python scripts directory
pydir=${PWD}/python

# Common data directory
datadir=${PWD}/data

# Prefix
prefix=full

# Geometry
if test ${prefix} = "full"; then
  nx=1796
  ny=2321
  lonCentre=15.001164907304375
  latCentre=63.50071790568658
  lon0=15.0
  lat0=63.0
elif test ${prefix} = "small"; then
  nx=200
  ny=200
  lonCentre=9.5
  latCentre=56
  lon0=9.5
  lat0=56
else
  echo "Wrong domain prefix, exiting"
fi

# Standard-deviations
sigmab=0.39 #0.316
sigmao=0.35 #0.1

# Screening factor
scrFac=3.0

# Correlation length-scales
rh=30e3
rv=1500.0

# Create directory
mkdir -p ${prefix}

# Clean log
rm -f ${prefix}/log

# Get data
#python ${pydir}/getData.py ${datadir}/fine_grid.nc ${datadir}/background.nc ${datadir}/observations_all.nc

# Restrict background to the domain
update_json restrictBackground.json
mpirun -n ${np} ${bindir}/quench_convertstate.x ${prefix}/restrictBackground.json

# Restrict observations to the domain
python ${pydir}/restrictObservations.py ${prefix}/background.nc ${datadir}/observations_all.nc ${prefix}/observations_all.nc

# Control observations quality
python ${pydir}/controlObservations.py ${prefix}/background.nc ${prefix}/observations_all.nc ${prefix}/observations.nc ${sigmab} ${sigmao} ${scrFac}

# Hollingworgh-Lönnberg diagnostics
#python ${pydir}/diagnoseHL.py ${prefix}/observations.nc

# Prepare background error
# - Standard-deviation
update_json createStdDev.json
mpirun -n ${np} ${bindir}/quench_convertstate.x ${prefix}/createStdDev.json
# - BUMP
#update_json prepareBackgroundErrorBUMP.json
#mpirun -n ${np} ${bindir}/quench_error_covariance_toolbox.x ${prefix}/prepareBackgroundErrorBUMP.json
#update_json diracBackgroundErrorBUMP.json
#mpirun -n ${np} ${bindir}/quench_error_covariance_toolbox.x ${prefix}/diracBackgroundErrorBUMP.json
# - FastLAM
update_json prepareBackgroundErrorFastLAM.json
mpirun -n ${np} ${bindir}/quench_error_covariance_toolbox.x ${prefix}/prepareBackgroundErrorFastLAM.json
update_json diracBackgroundErrorFastLAM.json
mpirun -n ${np} ${bindir}/quench_error_covariance_toolbox.x ${prefix}/diracBackgroundErrorFastLAM.json

for it in $(seq 1 ${nt}); do
  # Split observation
  python ${pydir}/splitObservations.py ${prefix}/observations.nc ${prefix}/observations_control.nc ${prefix}/observations_assim.nc

  # Run OI
  python ${pydir}/runOI.py ${prefix}/background.nc ${prefix}/observations_assim.nc ${prefix}/analysis_OI.nc ${prefix}/observations_assim_OI.nc
  ncdiff -O ${prefix}/analysis_OI.nc ${prefix}/background.nc ${prefix}/increment_OI.nc

  # Run 3DVar
  update_json 3dvarFastLAM.json
  mpirun -n ${np} ${bindir}/quench_variational.x ${prefix}/3dvarFastLAM.json
  ncdiff -O ${prefix}/analysis_3dvar.nc ${prefix}/background.nc ${prefix}/increment_3dvar.nc

  # Run Desroziers' diagnostics
#  python ${pydir}/diagnoseDesroziers.py ${prefix}/observations_assim_3dvar.nc

  # Compute errors with control observations
  python ${pydir}/computeErrors.py ${prefix}/observations_control.nc ${prefix}/background.nc ${prefix}/analysis_OI.nc ${prefix}/analysis_3dvar.nc ${prefix}/errors_${it}.nc
done

# Compare OI/3dvar
for it in $(seq 2 ${nt}); do
  ncrcat -A ${prefix}/errors_${it}.nc ${prefix}/errors.nc
done
python ${pydir}/gatherErrors.py ${prefix}/errors.nc ${prefix}/background.nc
