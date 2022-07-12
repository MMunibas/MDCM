#!/bin/bash

# This script finds the best combination of fragment results that combine to a
# given number of molecular charges. A total RMSE is estimated by averaging 
# over the individual fragment RMSEs, the combination of fragment models with
# lowest predicted overall RMSE is returned by this script for each given
# number of charges

# FOLDERS
ROOTDIR=/home/devereux/MDCM-release/examples/fragment-multiconf
WORKDIR=$ROOTDIR/4-combine-frags
BINDIR=$ROOTDIR/../../bin
FRAGDIR=$ROOTDIR/3-fit-frags # folder containing fitted fragments
# FITTING PARAMETERS
NFRAG=2             #number of fragments fitted
NFIT=5              #number of fits per fragment
MINCHGS=18          #minimum number of charges for whole molecule
MAXCHGS=36          #maximum number of charges for whole molecule
MINFRAGCHGS=8       #minimum charges used to fit fragments in previous step
MAXFRAGCHGS=20      #maximum charges used to fit fragments in previous step


python3 $BINDIR/find-min-rmse.py $BINDIR $FRAGDIR $WORKDIR $NFRAG $NFIT $MINCHGS $MAXCHGS $MINFRAGCHGS $MAXFRAGCHGS


