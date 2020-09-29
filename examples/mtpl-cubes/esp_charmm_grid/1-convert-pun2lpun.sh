#!/bin/bash

# requires rdkit, installed via anaconda:
# $ conda create -c rdkit -n my-rdkit-env rdkit
# $ conda activate my-rdkit-env

# This script converts a GDMA punch file to CHARMM lpun format

WORKDIR=/home/devereux/eric/esp_charmm_grid
REFDIR=$WORKDIR/ref
BINDIR=$WORKDIR/bin
#INPUTS
XYZFILE=$REFDIR/acrolein_1200.xyz
PUNFILE=$REFDIR/acrolein_1200.pun
#OUTPUTS
SDFFILE=$REFDIR/acrolein_1200.sdf
LRAPUNFILE=$REFDIR/acrolein_1200_lra.pun

cd $WORKDIR

# generate an SDF format file from an xyz
# coordinates should match those in the GDMA punch file
babel -ixyz $XYZFILE -osdf > $SDFFILE

# activate conda rdkit environment
. ~/.bashrc
conda activate my-rdkit-env
# add local axes to punch file
python3 $BINDIR/calc_LRA.py -in $SDFFILE -punfile $PUNFILE -lpunfile $LRAPUNFILE -punxyz
# convert gdma punch file to CHARMM lpun file in local axis system
python3 $BINDIR/pun2charmmlpun.py -pun $LRAPUNFILE
OUTNAME=$(echo $LRAPUNFILE | sed "s/\.pun/\.lpun/g")
echo "lpun file $OUTNAME written"
