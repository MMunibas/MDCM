#!/bin/bash

# Here we fit new atomic multipoles to the reference MEP from the PCUBE cube
# file of each conformer separately. Atomic charges are constrained to be the
# same for all conformers.

# In the second step we generate a Gaussian-format MEP cube file for each
# conformer using our newly fitted multipoles

# In the final step we compare the MEP of each new cube file to the reference
# MEP from our original cube file at different distances from the nuclei, to
# evaluate the quality of the fit

# FOLDERS
ROOT=/home/devereux/MDCM-release/examples/fragment-multiconf
WORKDIR=$ROOT/1-mtp-fit
BINDIR=$ROOT/../../bin
CUBEDIR=$ROOT/ref
# INPUT
PCUBE=$CUBEDIR/equi.pot.mod.cube    # ESP CUBE file supplied by user
DCUBE=$CUBEDIR/equi.dens.mod.cube   # Density CUBE file supplied by user 
PCUBE2=$CUBEDIR/conf2.pot.mod.cube  # ESP CUBE file supplied by user (confomer 2)
DCUBE2=$CUBEDIR/conf2.dens.mod.cube # Density CUBE file supplied by user (conformer 2)
# OUTPUT
MTPFILE=$WORKDIR/fitted-mtpl-equi.dat
MTPFILE2=$WORKDIR/fitted-mtpl-conf2.dat

cd $WORKDIR

# CONFORMER 1

# Fit new multipoles to reference MEP from $PCUBE 
$BINDIR/mtpfit.py -pot $PCUBE -dens $DCUBE -lmax 5 -qtot 0.0
mv fitted-mtpl.dat $MTPFILE
# Generate cube file from the newly fitted multipoles
$BINDIR/pcubefit.x -v -generate -multipole -esp $PCUBE -dens $DCUBE -mtpfile $MTPFILE > generate-cube-l5-equi.log
# Examine quality of fitted multipoles by comparing newly fitted multipolar and reference
# MEP
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE -esp2 $WORKDIR/ditriantapole_expansion.cube -dens $DCUBE > analyze-cube-l5-equi.log

# CONFORMER 2

# Second conformer: for subsequent fragment fitting the total fragment charges must
# be the same for both conformers, so we fix the atomic charges
$BINDIR/mtpfit.py -pot $PCUBE2 -dens $DCUBE2 -lmax 5 -qtot 0.0 -fixq $MTPFILE
mv fitted-mtpl.dat $MTPFILE2
# Generate cube file from the newly fitted multipoles
$BINDIR/pcubefit.x -v -generate -multipole -esp $PCUBE2 -dens $DCUBE2 -mtpfile $MTPFILE2 > generate-cube-l5-conf2.log
# Examine quality of fitted multipoles by comparing newly fitted multipolar and reference
# MEP
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE2 -esp2 $WORKDIR/ditriantapole_expansion.cube -dens $DCUBE2 > analyze-cube-l5-conf2.log

