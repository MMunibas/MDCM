#!/bin/bash

# This step contains a sample CHARMM setup for the residue TRP. It is doubly
# patched with neutral capping groups to create a reference system suitable
# for subsequent ab initio calculations, with an ESP that is not dominated by 
# charged capping groups. Note that some minor corrections to duplicate atom
# types (if an atom type appears in both terminal groups) and to add missing
# parameters may be necessary to allow double-patching a single residue in
# CHARMM.
#
# The atomic charges are dumped by CHARMM and extracted below, so that the 
# charges of the backbone atoms and capping groups can be subtracted from the
# reference ESP and frozen during the sidechain fit.
#
# As well as the minimized geometry, a second, distored geometry was selected
# using VMD to provide a second conformer to fit the ESP using conformational
# averaging (here named "frame3.pdb").

REFCHARMM=/opt/cluster/programs/charmm/developer/dev-release-dcm/build/cmake/charmm
NPROC=4
JOB=charmm-patch

ulimit -s 10240

module load gcc/gcc-9.2.0-openmpi-4.0.2-ib

# run gas-phase simulation
mpirun -np $NPROC $REFCHARMM -i ${JOB}.inp -o ${JOB}.log

# extract charges from output file
grep "( LIG  TRP  1" ${JOB}.log  | awk '{print $5,$7}' > charges.dat

