#!/bin/bash

# In this step we test the newly fitted sidechain model with two CHARMM runs.
# The first simply dumps the coordinates of all nuclei and MDCM charges for the
# generated tripeptide. The positions of charges relative to nuclei should be 
# the same was obtained during fitting.
#
# The second input runs basic dynamics with the model as a stability check and
# for direct comparison with the original CHARMM sidechain model if desired.

#DIRECTORIES
ROOT=/home/devereux/MDCM-release/examples/sidechain-fit
WORKDIR=$ROOT/10-charmm-test

#PARAMETERS
CHARMM=/opt/cluster/programs/charmm/developer/dev-release-dcm/build/cmake/charmm
NPROC=4
JOB1=charmm-trp-dcm
JOB2=charmm-trp-dcm-dyna

ulimit -s 10240
cd $WORKDIR

module load gcc/gcc-9.2.0-openmpi-4.0.2-ib

# run test to dump generated coordinates and dcm charge positions
mpirun -np $NPROC $CHARMM -i ${JOB1}.inp -o ${JOB1}.log

# run gas-phase simulation
mpirun -np $NPROC $CHARMM -i ${JOB2}.inp -o ${JOB2}.log


