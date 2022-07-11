#!/bin/bash

# Here the ab initio reference data is generated using Gaussian 16. In a
# previous step (not shown here) the coordinates from the last step were used
# to run single-point energy evaluations with the checkpoint files specified 
# below. These files are used here to generate the CUBE files that are needed
# for ESP fitting

WORKDIR=/home/devereux/MDCM-release/examples/sidechain-fit/2-cube-files
GAUSSDIR=/opt/cluster/programs/g16-c.01/g16/
CHK1=equi.chk
CHK2=conf2.chk


cd $WORKDIR

# create formatted checkpoint files:
$GAUSSDIR/formchk $CHK1
$GAUSSDIR/formchk $CHK2

FCHK1=$(basename $CHK1 .chk).fchk
FCHK2=$(basename $CHK2 .chk).fchk
DCUBE1=$(basename $CHK1 .chk).dens.cube
DCUBE2=$(basename $CHK2 .chk).dens.cube
PCUBE1=$(basename $CHK1 .chk).pot.cube
PCUBE2=$(basename $CHK2 .chk).pot.cube

# Create CUBE files for 1st conformer (whole system)
$GAUSSDIR/cubegen 0 density $FCHK1 $DCUBE1 -2 h
$GAUSSDIR/cubegen 0 potential $FCHK1 $PCUBE1 -2 h

# Create CUBE files for 2nd conformer (whole system)
$GAUSSDIR/cubegen 0 density $FCHK2 $DCUBE2 -2 h
$GAUSSDIR/cubegen 0 potential $FCHK2 $PCUBE2 -2 h
