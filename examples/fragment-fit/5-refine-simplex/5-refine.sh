#!/bin/bash

# In this step we refine the best guess created in the previous step by combining
# the separate fragment fits. The script accepts a range of charge models to be
# refined ranging from $MINCHGS to $MAXCHGS for the whole molecule. For larger
# systems with more charges it is necessary to perform refinement using a simplex
# algorithm, as done here for demonstration purposes, or at least to use a
# relaxed convergence criterion with the -conv flag (see manual). For smaller
# systems an additional DE refinement for the whole molecule may yield better
# results.

# Choose the "PCUBE" reference MEP file as appropriate 

# Modify Slurm script below for your cluster environment

# FOLDERS
ROOT=/home/devereux/MDCM-release/examples/fragment-fit
WORKDIR=$ROOT/5-refine-simplex
BINDIR=$ROOT/../../bin
REFDIR=$ROOT/../ref
FRAGDIR=$ROOT/4-combine-frags  # folder containing combined fragment xyz files
# INPUT
MTPFILE=$ROOT/1-mtp-fit/fitted-mtpl.dat
PCUBE=$REFDIR/benzene-pot.cube   # ESP reference data
DCUBE=$REFDIR/benzene-dens.cube  
# FITTING PARAMETERS
MINNTRY=1      # tries for largest number of charges (DE only)
MAXNTRY=1      # tries for smallest number of charges (DE only)
MINCHGS=12     # minimum number of charges for whole molecule
MAXCHGS=24     # maximum number of charges for whole molecule
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

$BINDIR/pcubefit.x -xyz $INITIALXYZ -mtpfile $MTPFILE -simplex -esp $PCUBE -dens $DCUBE -nacmax $MAXATMCHG -ntry $NTRY -v > $i-charges.out
" > $REFINEDIR/$i-charges/$i-charges.sh
  sbatch $REFINEDIR/$i-charges/$i-charges.sh


done

