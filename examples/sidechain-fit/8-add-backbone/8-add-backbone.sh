#!/bin/bash

# In this step the frozen charges are added onto the fitted sidechain charges
# to allow evaluation of the performance of the full model, compared to the 
# reference ab initio CUBE files of the full system.

# DIRECTORIES
ROOT=/home/devereux/MDCM-release/examples/sidechain-fit
WORKDIR=$ROOT/8-add-backbone
CUBEDIR=$ROOT/2-cube-files
BINDIR=$ROOT/../../bin/

#INPUTS
CHGFILE=$ROOT/1-charmm-files/charges.dat # created in step 1
XYZFILE=$ROOT/6-fit-molecule/fit3/21charges.xyz # chosen fitted model
DCUBE=$CUBEDIR/equi.dens.cube # cube file needed for nuclear coorindates in same global axis as fitted charges
PCUBE=$CUBEDIR/equi.pot.cube # full ESP cube file (including frozen charges) for analysis
FRAMEFILE=$ROOT/6-fit-molecule/trp-frames.txt
ATOMLIST="1,2,3,4,5,6,7,26,27" # the list of atoms that were frozen during fitting

cd $WORKDIR

# This script adds the charges that were frozen during fitting to the newly fitted 
# MDCM of the sidechain
$BINDIR/add_charges_to_xyz.pl $XYZFILE $CHGFILE $DCUBE $ATOMLIST > combined.xyz

# generate a new CUBE file from the full system charge model
$BINDIR/pcubefit.x -generate -xyz combined.xyz -esp $PCUBE -dens $DCUBE -v

NAME=$(ls -1rt *charges.cube | sed 's/.cube//g')

# compare that CUBE file to the original ab initio reference
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE -esp2 *charges.cube -dens $DCUBE > analyze-cube-${NAME}.log

