#!/bin/bash

# This script finds the best combination of fragment results that combine to a
# given number of molecule charges. A total RMSE is estimated by averaging 
# over the individual fragment RMSEs, the combination of fragment models with
# lowest predicted overall RMSE is returned by this script for each given
# number of charges

# FOLDERS
ROOT=/home/devereux/MDCM-release/examples/fragment-fit
WORKDIR=$ROOT/4-combine-frags
BINDIR=/home/devereux/MDCM-release/bin
FRAGDIR=$ROOT/3-fit-fragments # folder containing fitted fragments
# FITTING PARAMETERS
NFRAG=2             #number of fragments fitted
NFIT=5              #number of fits per fragment
MINCHGS=12          #minimum number of charges for whole molecule
MAXCHGS=24          #maximum number of charges for whole molecule
MINFRAGCHGS=6       #minimum charges used to fit fragments in previous step
MAXFRAGCHGS=12      #maximum charges used to fit fragments in previous step


python3 $BINDIR/find-min-rmse.py $BINDIR $FRAGDIR $WORKDIR $NFRAG $NFIT $MINCHGS $MAXCHGS $MINFRAGCHGS $MAXFRAGCHGS


