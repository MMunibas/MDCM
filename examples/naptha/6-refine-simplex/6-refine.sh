#!/bin/bash

# In this step we refine the best guess created in the previous step by combining
# the separate fragment fits. The script accepts a range of charge models to be
# refined ranging from $MINCHGS to $MAXCHGS for the whole molecule. For larger
# systems with more than ca. 18 charges for the whole molecule it is necessary to
# perform refinement using a simplex algorithm, as done here. For smaller systems
# an additional DE refinement for the whole molecule may yield better results.

# Choose the "PCUBE" reference MEP file as appropriate (fit to ab initio reference,
# multipolar MEP reference etc.)

# Modify Slurm script below for your cluster environment

# FOLDERS
ROOT=/home/devereux/naphta04/mdcm_fit_to_mtp/
WORKDIR=$ROOT/6-refine-simplex
BINDIR=$ROOT/bin
REFDIR=$WORKDIR/ref
FRAGDIR=$ROOT/5-combine-frags  # folder containing combined fragment xyz files
# INPUT
MTPFILE=$ROOT/2-mtp-fit/fitted-mtpl.dat
PCUBE=$ROOT/2-mtp-fit/ditriantapole_expansion.cube # used as reference data for the fit
DCUBE=$REFDIR/naphta04.dens.cube
# FITTING PARAMETERS
MINNTRY=1      #tries for largest number of charges (DE only)
MAXNTRY=1      #tries for smallest number of charges (DE only)
MINCHGS=30     #minimum number of charges for whole molecule
MAXCHGS=64     #maximum number of charges for whole molecule
MAXATMCHG=4    # max charges per atom used in atom fits (step 3)

REFINEDIR=$WORKDIR/refined

mkdir -p $REFINEDIR
cd $REFINEDIR

qrange=$(echo "$MAXCHGS - $MINCHGS" | bc -l)
tryrange=$(echo "$MAXNTRY - $MINNTRY" | bc -l)

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

$BINDIR/pcubefit.x -xyz $INITIALXYZ $MTPFILE -simplex -esp $PCUBE -dens $DCUBE -nacmax $MAXATMCHG -ntry $NTRY -v > $i-charges.out
" > $REFINEDIR/$i-charges/$i-charges.sh
  sbatch $REFINEDIR/$i-charges/$i-charges.sh


done
