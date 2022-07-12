#!/bin/bash

# This script evaluates the quality of a fitted model specified by the
# user as "XYZFILE". The model is used to generate an ESP CUBE file for a
# given conformer, which is then compared to a reference CUBE file for that
# conformer, for example from an ab initio computation.
#
# As we have more than one conformer, we need to define local reference 
# axes to describe how charges are displaced with conformational change

# FOLDERS
ROOTDIR=/home/devereux/MDCM-release/examples/fragment-multiconf
WORKDIR=$ROOTDIR/6-analysis
BINDIR=$ROOTDIR/../../bin     # path to MDCM scripts and binaries
REFDIR=$ROOTDIR/ref
# INPUT
XYZFILE=$ROOTDIR/5-refine/refined/26-charges/26_charges_refined.xyz # chosen fitted model
PCUBE=$REFDIR/equi.pot.mod.cube     # reference ESP data
DCUBE=$REFDIR/equi.dens.mod.cube    # reference electron density
PCUBE2=$REFDIR/conf2.pot.mod.cube   # reference ESP data (conformer 2)
DCUBE2=$REFDIR/conf2.dens.mod.cube  # reference electron density (conformer 2)
FRAMEFILE=$REFDIR/trp-frames.txt    # reference axis file supplied by user

cd $WORKDIR

NAME=$(basename $XYZFILE) # strips path
NAME=${NAME%%_*} # strips _charges_refined.xyz from filename
NAME=${NAME}charges # adds "charges" to end of name

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

