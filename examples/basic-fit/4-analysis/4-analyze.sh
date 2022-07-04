#!/bin/bash

# This script evaluates the quality of a fitted model specified by the
# user as "XYZFILE". The model is used to generate an ESP CUBE file, which
# is then compared to a reference CUBE file, for example from an ab initio
# computation.

# FOLDERS
ROOTDIR=/home/devereux/MDCM-release/examples/basic-fit
WORKDIR=$ROOTDIR/4-analysis
BINDIR=/home/devereux/MDCM-release/bin
REFDIR=$ROOTDIR/../ref
# INPUT
XYZFILE=$ROOTDIR/3-fit-molecule/fit3/7charges.xyz # chosen fitted model
PCUBE=$REFDIR/h2o-pot.cube  # reference ESP CUBE file
DCUBE=$REFDIR/h2o-dens.cube # reference density CUBE file

cd $WORKDIR

NAME=$(basename $XYZFILE)
NAME=${NAME%.*}

# generate a new CUBE file from the charge model
$BINDIR/pcubefit.x -generate -xyz $XYZFILE -esp $PCUBE -dens $DCUBE -v

# compare that CUBE file to a reference
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE -esp2 ${NAME}.cube -dens $DCUBE > analyze-cube-${NAME}.log

