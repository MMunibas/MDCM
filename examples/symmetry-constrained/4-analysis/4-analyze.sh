#!/bin/bash

# In this step we create ESP cube files from the fitted charge models to compare with
# the reference ESP 

# FOLDERS
ROOTDIR=/home/devereux/MDCM-release/examples/symmetry-constrained
WORKDIR=$ROOTDIR/4-analysis
BINDIR=/home/devereux/MDCM-release/bin
REFDIR=$ROOTDIR/../ref
FITDIR=$ROOTDIR/3-fit-molecule
# INPUT
XYZFILE=$FITDIR/fit1/30charges.xyz # preferred fitted model
PCUBE=$REFDIR/benzene-pot.cube     # reference ESP CUBE file
DCUBE=$REFDIR/benzene-dens.cube    # reference density CUBE file

cd $WORKDIR

NAME=$(basename $XYZFILE)
NAME=${NAME%.*}

$BINDIR/pcubefit.x -generate -xyz $XYZFILE -esp $PCUBE -dens $DCUBE -v > ${NAME}.out

# Examine quality of fitted charges by comparing newly fitted model and reference
# MEP
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE -esp2 ${NAME}.cube -dens $DCUBE > analyze-cube-${NAME}.log

