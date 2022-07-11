#!/bin/bash

# Here we fit new atomic multipoles to the reference MEP from the PCUBE cube file
# of the equilibrium geometry only (this CUBE file has had the backbone atoms and
# capping groups removed in a previous step)

# In the second step we generate a Gaussian-format MEP cube file using our newly
# fitted multipoles

# In the final step we compare the MEP in the new cube file to the reference MEP
# from our original cube file at different distances from the nuclei, to evaluate
# the quality of the fit

# FOLDERS
ROOT=/home/devereux/MDCM-release/examples/sidechain-fit
WORKDIR=$ROOT/4-mtp-fit
BINDIR=$ROOT/../../bin
CUBEDIR=$ROOT/3-modify-cubes
# INPUT
PCUBE=$CUBEDIR/equi.pot.mod.cube  # ESP CUBE file supplied by user
DCUBE=$CUBEDIR/equi.dens.mod.cube # Density CUBE file supplied by user 
# OUTPUT
MTPFILE=$WORKDIR/fitted-mtpl.dat

cd $WORKDIR

# Fit new multipoles to reference MEP from $PCUBE 
$BINDIR/mtpfit.py -pot $PCUBE -dens $DCUBE -lmax 5 -qtot 0.0

# Generate cube file from the newly fitted multipoles
$BINDIR/pcubefit.x -v -generate -multipole -esp $PCUBE -dens $DCUBE -mtpfile $MTPFILE > generate-cube-l5.log

# Examine quality of fitted multipoles by comparing newly fitted multipolar and reference
# MEP
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE -esp2 $WORKDIR/ditriantapole_expansion.cube -dens $DCUBE > analyze-cube-l5.log

