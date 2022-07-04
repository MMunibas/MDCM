#!/bin/bash

# Here we fit new atomic multipoles to the reference MEP from the PCUBE cube file

# In the second step we generate a Gaussian-format MEP cube file using our newly
# fitted multipoles

# In the final step we compare the MEP in the new cube file to the reference MEP
# from our original cube file at different distances from the nuclei, to evaluate
# the quality of the fit

# FOLDERS
ROOT=/home/devereux/MDCM-release/examples/basic-fit
WORKDIR=$ROOT/1-mtp-fit
BINDIR=$ROOT/../../bin
REFDIR=$ROOT/../ref
# INPUT
PCUBE=$REFDIR/h2o-pot.cube  # ESP CUBE file supplied by user
DCUBE=$REFDIR/h2o-dens.cube # Density CUBE file supplied by user 
# OUTPUT
MTPFILE=$WORKDIR/fitted-mtpl.dat

cd $WORKDIR

# Fit new multipoles to reference MEP from $PCUBE (can be ab initio or multipolar MEP)
$BINDIR/mtpfit.py -pot $PCUBE -dens $DCUBE -lmax 5 -qtot 0.0

# Generate cube file from the newly fitted multipoles
$BINDIR/pcubefit.x -v -generate -multipole -esp $PCUBE -dens $DCUBE -mtpfile $MTPFILE > generate-cube-l5.log

# Examine quality of fitted multipoles by comparing newly fitted multipolar and reference
# MEP
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE -esp2 $WORKDIR/ditriantapole_expansion.cube -dens $DCUBE > analyze-cube-l5.log

