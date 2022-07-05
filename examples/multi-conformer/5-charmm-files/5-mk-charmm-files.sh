#!/bin/bash

# This script converts a fitted charge model to a ".dcm" file that can be used
# in CHARMM. Pay attention to the atom ordering in the "frames.txt" file.
# assuming that the atom ordering used for fitting matches the order in the 
# CHARMM topology file, the frames.txt file contains groups of 3 atoms that make
# up each frame. The label "BO" or "BI" requests either a bond as the z-axis or
# the A-B-C bisector angle as the z-axis for the given frame.
#
# Example:
# SCN                        !Title
# 2    1    3    BO          !Frame is defined by atoms 2,1,3 using bonds for z-axes
#

# FOLDERS
ROOTDIR=/home/devereux/MDCM-release/examples/multi-conformer
WORKDIR=$ROOTDIR/5-charmm-files
BINDIR=/home/devereux/MDCM-release/bin
REFDIR=$ROOTDIR/../ref
FITDIR=$ROOTDIR/3-fit-molecule
# INPUT
XYZFILE=$FITDIR/fit3/7charges.xyz # fitted model
DCUBE=$REFDIR/h2o-dens.cube
FRAMEFILE=$REFDIR/h2o-frames.txt

cd $WORKDIR

NAME=$(basename $XYZFILE)
NAME=${NAME%.*}

$BINDIR/comb-xyz-to-dcm.pl $XYZFILE $DCUBE $FRAMEFILE ${NAME}.dcm


