#!/bin/bash

# In this step the reference cube files generated using Gaussian are processed
# by a script to remove the ESP contribution and the associated atoms listed in
# the array "ATOMLIST". These atoms correspond to the backbone and capping
# group atoms that will be frozen during sidechain fitting.

# DIRECTORIES
ROOT=/home/devereux/MDCM-release/examples/sidechain-fit
CUBEDIR=$ROOT/2-cube-files      # constains reference CUBE files of whole system
WORKDIR=$ROOT/3-modify-cubes    
BINDIR=$ROOT/../../bin/         # contains MDCM scripts and code
# INPUT
CHGFILE=$ROOT/1-charmm-files/charges.dat
DCUBE1=$CUBEDIR/equi.dens.cube  # reference electron density of full system
PCUBE1=$CUBEDIR/equi.pot.cube   # reference ESP of full system
DCUBE2=$CUBEDIR/conf2.dens.cube # reference density of second conformer
PCUBE2=$CUBEDIR/conf2.pot.cube  # reference ESP of second conformer
ATOMLIST="1,2,3,4,5,6,7,26,27"  # list of atom indices to freeze during fitting

cd $WORKDIR

# 1st conformer: remove ESP due to frozen charges and write new CUBE file
$BINDIR/remove_frag_esp_from_cube.pl $CHGFILE $DCUBE1 $PCUBE1 $ATOMLIST
PCUBE=$(basename $PCUBE1 .cube).mod.cube
DCUBE=$(basename $DCUBE1 .cube).mod.cube
mv modified.dens.cube $DCUBE
mv modified.pot.cube $PCUBE

# 2nd conformer: remove ESP due to frozen charges and write new CUBE file
$BINDIR/remove_frag_esp_from_cube.pl $CHGFILE $DCUBE2 $PCUBE2 $ATOMLIST
PCUBE=$(basename $PCUBE2 .cube).mod.cube
DCUBE=$(basename $DCUBE2 .cube).mod.cube
mv modified.dens.cube $DCUBE
mv modified.pot.cube $PCUBE

