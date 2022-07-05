#!/bin/bash

# This script evaluates the quality of a fitted model specified by the
# user as "XYZFILE". The model is used to generate an ESP CUBE file for a
# given conformer, which is then compared to a reference CUBE file for that
# conformer, for example from an ab initio computation.
#
# As we have more than one conformer, we need to define local reference 
# axes to describe how charges are displaced with conformational change

# FOLDERS
ROOTDIR=/home/devereux/MDCM-release/examples/multi-conformer
WORKDIR=$ROOTDIR/4-analysis
BINDIR=/home/devereux/MDCM-release/bin
REFDIR=$ROOTDIR/../ref
# INPUT
XYZFILE=$ROOTDIR/3-fit-molecule/fit3/7charges.xyz # chosen fitted model
PCUBE=$REFDIR/h2o-pot.cube  # reference ESP CUBE file
DCUBE=$REFDIR/h2o-dens.cube # reference density CUBE file
PCUBE2=$REFDIR/h2o-pot-acute.cube  # reference ESP CUBE file (conformer 2)
DCUBE2=$REFDIR/h2o-dens-acute.cube # reference density CUBE file (conformer 2)
PCUBE3=$REFDIR/h2o-pot-obtuse.cube  # reference ESP CUBE file (conformer 3)
DCUBE3=$REFDIR/h2o-dens-obtuse.cube # reference density CUBE file (conformer 3)
FRAMES=$REFDIR/h2o-frames.txt

cd $WORKDIR

NAME=$(basename $XYZFILE)
NAME=${NAME%.*}

# First assess quality with respect to conformer 1, which has the same global
# frame as the fitted MDCM in XYZFILE so doesn't need to be transformed

# generate a new CUBE file from the charge model
$BINDIR/pcubefit.x -generate -xyz $XYZFILE -esp $PCUBE -dens $DCUBE -v

# compare that CUBE file to a reference
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE -esp2 ${NAME}.cube -dens $DCUBE > analyze-cube-${NAME}.log

# Now so the same for conformer 2, this time transforming the charge model 
# using local axes to follow the new molecular conformation

# generate a new CUBE file from the charge model
$BINDIR/pcubefit.x -generate -xyz $XYZFILE -esp $PCUBE -dens $DCUBE -esp $PCUBE2 -dens $DCUBE2 -frames $FRAMES -v

# compare that CUBE file to a reference
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE2 -esp2 ${NAME}.cube -dens $DCUBE2 > analyze-cube-${NAME}-2.log

# Now so the same for conformer 3, again transforming the charge model 
# using local axes to follow the new molecular conformation

# generate a new CUBE file from the charge model
$BINDIR/pcubefit.x -generate -xyz $XYZFILE -esp $PCUBE -dens $DCUBE -esp $PCUBE3 -dens $DCUBE3 -frames $FRAMES -v

# compare that CUBE file to a reference
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE3 -esp2 ${NAME}.cube -dens $DCUBE3 > analyze-cube-${NAME}-3.log


