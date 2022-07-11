#!/bin/bash

# This script evaluates the quality of a fitted model specified by the
# user as "XYZFILE". The model is used to generate an ESP CUBE file for a
# given conformer, which is then compared to a reference CUBE file for that
# conformer, for example from an ab initio computation.
#
# As we have more than one conformer, we need to define local reference 
# axes to describe how charges are displaced with conformational change. In
# this example the necessary "frames" file was already created for the previous
# step.
#
# Note that we still work with the truncated system, where the frozen atoms and
# associated ESP have been removed from the reference CUBE files and hence the
# quality of the fit refers to the sidechain only.

# FOLDERS
ROOTDIR=/home/devereux/MDCM-release/examples/sidechain-fit
WORKDIR=$ROOTDIR/7-analysis
BINDIR=$ROOTDIR/../../bin
CUBEDIR=$ROOTDIR/3-modify-cubes/
# INPUT
XYZFILE=$ROOTDIR/6-fit-molecule/fit3/21charges.xyz # chosen fitted model
PCUBE=$CUBEDIR/equi.pot.mod.cube    
DCUBE=$CUBEDIR/equi.dens.mod.cube
PCUBE2=$CUBEDIR/conf2.pot.mod.cube 
DCUBE2=$CUBEDIR/conf2.dens.mod.cube
FRAMEFILE=$ROOTDIR/6-fit-molecule/trp-frames.txt

cd $WORKDIR

NAME=$(basename $XYZFILE) # strips path
NAME=${NAME%.*} # strips .xyz from filename

# First assess quality with respect to conformer 1, which has the same global
# frame as the fitted MDCM in XYZFILE so doesn't need to be transformed

# generate a new CUBE file from the charge model
$BINDIR/pcubefit.x -generate -xyz $XYZFILE -esp $PCUBE -dens $DCUBE -v

# compare that CUBE file to a reference
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE -esp2 ${NAME}.cube -dens $DCUBE > analyze-cube-${NAME}.log

# Now so the same for conformer 2, this time transforming the charge model 
# using local axes to follow the new molecular conformation

# generate a new CUBE file from the charge model
$BINDIR/pcubefit.x -generate -xyz $XYZFILE -esp $PCUBE -dens $DCUBE -esp $PCUBE2 -dens $DCUBE2 -frames $FRAMEFILE -v

# compare that CUBE file to a reference
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE2 -esp2 ${NAME}.cube -dens $DCUBE2 > analyze-cube-${NAME}-2.log

