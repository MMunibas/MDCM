#!/bin/bash

# In this step we create ESP cube files from the fitted charge models to compare with
# the reference ESP visually

# FOLDERS
ROOTDIR=/pchem-data/meuwly/devereux/MDCM/examples/benzene-sym
WORKDIR=$ROOTDIR/5-analysis
BINDIR=/pchem-data/meuwly/devereux/MDCM/bin
REFDIR=$ROOTDIR/ref
FITDIR=$ROOTDIR/4-fit-molecule
# INPUT
XYZFILE=$FITDIR/fit1/24charges.xyz # fitted model
PCUBE=$REFDIR/benzene-pot.cube # again, not really used for fragments
DCUBE=$REFDIR/benzene-dens.cube

cd $WORKDIR
module load gcc/gcc-7.4.0-openmpi-3.1.4-cuda-for-11.4

NAME=$(basename $XYZFILE)
NAME=${NAME%.*}

$BINDIR/cudacubefit.x -generate -xyz $XYZFILE -esp $PCUBE -dens $DCUBE -v > ${NAME}.out

# Examine quality of fitted charges by comparing newly fitted model and reference
# MEP
$BINDIR/cudacubefit.x -v -analysis -esp $PCUBE -esp2 ${NAME}.cube -dens $DCUBE > analyze-cube-${NAME}.log

