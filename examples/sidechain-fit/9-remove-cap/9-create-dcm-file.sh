#!/bin/bash

# In this step we need to modify the full system, which includes capping groups
# and may have backbone charges that were modified during patching by CHARMM, 
# in order to remove the capping groups, reset the backbone charges to
# unpatched values and create a .dcm parameter file for the unpatched residue
# that can be used generically in CHARMM.
#
# This task is too complex to easily automate, the code below assumes that you 
# already modified the "atom_coords.xyz" file from the previous step by 
# removing the atoms that belong to the terminal groups and correcting the atom
# order to match the topology file. It also assumes that you removed the 
# charges from the "combined.xyz" MDCM file from the previous step to remove 
# the charges that belong to these atoms (the order of the charges isn't 
# important) and that you corrected any remaining backbone charges that were
# altered during patching (again using the CHARMM topology file as reference)
#
# Finally, a new frames file needs to be created for the unpatched residue,
# with atom indices corresponding to the topology file and updated atom_coords
# file

#DIRECTORIES
ROOT=/home/devereux/MDCM-release/examples/sidechain-fit
WORKDIR=$ROOT/9-remove-cap
BINDIR=$ROOT/../../bin

#INPUTS
CHGSFILE=$WORKDIR/combined_no_cap.xyz     # modified as described above
CRDSFILE=$WORKDIR/atom_coords_no_cap.xyz  # modified as described above
FRAMEFILE=$WORKDIR/trp-frames-no-cap.txt  # axis frames for unpatched residue (specified by user)

#OUTPUT
DCMFILE=$WORKDIR/trp.dcm

# Create .dcm parameter file for use in CHARMM
$BINDIR/comb-xyz-to-dcm-2.pl $CHGSFILE $CRDSFILE $FRAMEFILE $DCMFILE
