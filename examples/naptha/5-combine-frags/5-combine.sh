#!/bin/bash

# Script to find the best combination of fragment results that combine to a given number
# of molecule charges


# FOLDERS
ROOT=/home/devereux/MDCM-git/examples/naptha
WORKDIR=$ROOT/5-combine-frags
BINDIR=/home/devereux/MDCM-git/bin
FRAGDIR=$ROOT/4-fit-frags # folder containing fitted fragments
# FITTING PARAMETERS
NFRAG=5             #number of fragments fitted
NFIT=1              #number of fits per fragment
MINCHGS=18          #minimum number of charges for whole molecule
MAXCHGS=64          #maximum number of charges for whole molecule
MINFRAGCHGS=6       #minimum charges used to fit fragments in previous step
MAXFRAGCHGS=14      #maximum charges used to fit fragments in previous step


$BINDIR/find-min-rmse.py $ROOT $BINDIR $FRAGDIR $WORKDIR $NFRAG $NFIT $MINCHGS $MAXCHGS $MINFRAGCHGS $MAXFRAGCHGS
