#!/bin/bash

# This script runs a CHARMM input for a molecule interacting with a unit point
# charge, yielding the MEP across a grid

ROOTDIR=/home/devereux/eric/esp_charmm_grid
CHARMMDIR=$ROOTDIR/charmm
BINDIR=$ROOTDIR/bin
REFDIR=$ROOTDIR/ref
# make sure you set the pdb file in the CHARMM input to have the same coordinates as
# the reference cube, and set the cube grid variables in the input to match also
PCUBE=$REFDIR/acrolein_1200.chk.fchk.potential.cube
DCUBE=$REFDIR/acrolein_1200.chk.fchk.density.cube

module load charmm/c44b1-gcc4.8.5-openmpi3.0.0
cd $CHARMMDIR

# make sure you modify the CHARMM input file and supply the necessary coordinate
# files etc. defined therein before running this script. Note that the pdb file
# needs to have the "K" dummy atom appended to the bottom to evaluate the 
# interaction with a unit charge
echo "you modified ESP-template.inp to match $PCUBE already, right?"
echo
echo "the pdb file should have the same coords as $PCUBE (in Angstrom)"
echo
echo "the cube grid parameters in the CHARMM input should match $PCUBE"
echo
/opt/charmm/c44b1-gcc4.8.5 -i ESP-template.inp > /dev/null

# run script to generate a cube format file from the CHARMM output:
$BINDIR/mk_cube_from_charmm_dat.pl $PCUBE mep.dat

# finally, compare the CHARMM cube file to the reference cube:
$BINDIR/pcubefit.x -v -analysis -esp $PCUBE -esp2 $CHARMMDIR/charmm.cube -dens $DCUBE

