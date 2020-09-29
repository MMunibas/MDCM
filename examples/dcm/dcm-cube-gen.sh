#!/bin/bash

# This script converts a GDMA multipole file (maximum rank l=2, quadrupole) to
# an octahedral DCM representation with 6 charges per atom. As an accuracy check
# the DCM model is used to generate an MEP cube file which is compared to a
# reference GDMA multipolar cube file. The two should match closely.
# To obtain the GDMA cube file, see the first steps of the "naphta" example

# Note: the name of the file is determined by the residue name in the frames
# input file, i.e. $DCMFILE should match the residue name in $REFFRAMES

# Note: coordinates and atom numbering must be consistent for the cube files,
# GDMA punch file and 

#FOLDERS
ROOT=/home/devereux/MDCM-git/examples/dcm
REFDIR=$ROOT/ref
BINDIR=/home/devereux/MDCM-git/bin
#INPUT
PCUBE=$REFDIR/gdma-l2.cube            # Gaussian MEP cube file
DCUBE=$REFDIR/naphta.dens.cube        # Gaussian electron density cube file
GDMAFILE=$REFDIR/gdma-l2.punch        # GDMA multipoles implemented up to l=2 (quadrupole)
REFFRAMES=$REFDIR/naphta-frames.txt   # DCM reference axis systems for molecule
#OUPUT
DCMFILE=$ROOT/NAPHTA.dcm           # CHARMM format DCM parameter file

cd $ROOT

# create a .dcm parameter file from a reference GDMA punch file
# the name of the file is determined by the residue name in the frames input file
$BINDIR/punch-to-dcm.x $GDMAFILE $REFFRAMES

# generate a cube file from the DCM charge arrangement created in the 1st step
$BINDIR/mk_cube_from_dcm.pl $PCUBE $DCMFILE

# compare the cube file from DCM charges with the reference cube file
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE -esp2 $ROOT/dcm.cube -dens $DCUBE


