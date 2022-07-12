#!/bin/bash

# In this step we refine the best guess created in the previous step that combined
# the separate fragment fits. The script accepts a range of charge models to be
# refined ranging from $MINCHGS to $MAXCHGS for the whole molecule. For larger
# systems with more charges it is necessary to perform refinement using a simplex
# algorithm, as done here for demonstration purposes, or at least to use a
# relaxed convergence criterion with the -conv flag (see manual). For smaller
# systems an additional DE refinement for the whole molecule may yield slightly 
# better results.
#
# The simplex refinement uses both conformers as reference data.

# Modify Slurm script below for your cluster environment

# FOLDERS
ROOTDIR=/home/devereux/MDCM-release/examples/fragment-multiconf
WORKDIR=$ROOTDIR/5-refine
BINDIR=$ROOTDIR/../../bin
REFDIR=$ROOTDIR/ref
FRAGDIR=$ROOTDIR/4-combine-frags  # folder containing combined fragment xyz files
# INPUT
MTPFILE=$ROOTDIR/1-mtp-fit/fitted-mtpl-equi.dat # from step 1
PCUBE=$REFDIR/equi.pot.mod.cube    # reference ESP of whole system
DCUBE=$REFDIR/equi.dens.mod.cube   # reference electron density
MTPFILE2=$ROOTDIR/1-mtp-fit/fitted-mtpl-conf2.dat # from step 1
PCUBE2=$REFDIR/conf2.pot.mod.cube  # reference ESP of whole system (conformer 2)
DCUBE2=$REFDIR/conf2.dens.mod.cube # reference electron density (conformer 2)
FRAMEFILE=$REFDIR/trp-frames.txt   # multi-conformer fits require local reference axes
                                   # to be defined by the user
# FITTING PARAMETERS
MINNTRY=1      # tries for largest number of charges (DE only)
MAXNTRY=1      # tries for smallest number of charges (DE only)
MINCHGS=18     # minimum number of charges for whole molecule
MAXCHGS=36     # maximum number of charges for whole molecule
MAXATMCHG=4    # max charges per atom used in atom fits (step 3)

REFINEDIR=$WORKDIR/refined

mkdir -p $REFINEDIR
cd $REFINEDIR

qrange=$(echo "$MAXCHGS - $MINCHGS" | bc -l)    # range of charges 
# computational cost for a large number of charges is high, so we can 
# reduce the number of tries per fit for the largest models (DE only):
tryrange=$(echo "$MAXNTRY - $MINNTRY" | bc -l) 

# loop over all fits (different numbers of charges):
for ((i=$MINCHGS; i<=$MAXCHGS; i++)); do
  if [ ! -e $FRAGDIR/$i-combined.xyz ]; then
    echo "Error: file $FRAGDIR/$i-combined.xyz not found!"
    exit
  else
    INITIALXYZ="$FRAGDIR/$i-combined.xyz"
  fi
  # calculate NTRY for this molecule (linear decrease with increasing number of charges)
  if (( $MAXCHGS == $MINCHGS )); then  #avoid divide by zero
    NTRY=$MAXNTRY
  else
    tmp=$(echo "$i - $MINCHGS" | bc -l)
    tmp=$(echo "$MAXNTRY - $tryrange * $tmp / $qrange" | bc -l)
    NTRY=$(printf "%0.f" $tmp)
  fi
  mkdir -p $REFINEDIR/$i-charges
  echo "#!/bin/bash

#SBATCH --job-name=refine-${i}chg
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=long

cd $REFINEDIR/$i-charges

$BINDIR/pcubefit.x -xyz $INITIALXYZ -simplex -mtpfile $MTPFILE -esp $PCUBE -dens $DCUBE -mtpfile $MTPFILE2 -esp $PCUBE2 -dens $DCUBE2 -frames $FRAMEFILE -nacmax $MAXATMCHG -ntry $NTRY -v > $i-charges.out
" > $REFINEDIR/$i-charges/$i-charges.sh
  sbatch $REFINEDIR/$i-charges/$i-charges.sh


done

