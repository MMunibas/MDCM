#!/bin/bash

# Here we fit atomic multipoles to the reference MEP in the PCUBE cube file. No
# symmetry constraints are applied to the multipoles

# In the second step we generate a Gaussian-format MEP cube file using our newly
# fitted multipoles

# In the final step we compare the MEP in the new cube file to the reference MEP
# from our original cube file at different distances from the nuclei

# FOLDERS
ROOT=/home/devereux/MDCM-release/examples/symmetry-constrained
WORKDIR=$ROOT/1-mtp-fit
BINDIR=/home/devereux/MDCM-release/bin
REFDIR=$ROOT/../ref
# INPUT
PCUBE=$REFDIR/benzene-pot.cube # Reference ESP CUBE file supplied by user
DCUBE=$REFDIR/benzene-dens.cube # Reference electron density file supplied by user
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

