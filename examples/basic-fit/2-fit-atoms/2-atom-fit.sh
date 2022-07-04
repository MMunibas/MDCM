#!/bin/bash

# Here we fit atomic MDCM models to the atomic multipoles from the previous step.
# We fit models with 1,2,3..MAXCHG charges per atom, and for each fit we repeat
# NTRY times to increase chances of a good solution. We additionally run NFIT
# completely independent fits to more widely explore potential solutions. 

# Note: adapt the Slurm script below to whatever cluster you're using...

# FOLDERS
ROOT=/home/devereux/MDCM-release/examples/basic-fit
WORKDIR=$ROOT/2-fit-atoms
BINDIR=$ROOT/../../bin
REFDIR=$ROOT/../ref
MTPDIR=$ROOT/1-mtp-fit
# INPUT
DCUBE=$REFDIR/h2o-dens.cube
PCUBE=$REFDIR/h2o-pot.cube
MTPFILE=$MTPDIR/fitted-mtpl.dat # fitted multipoles from step 1
# FITTING PARAMETERS
NFIT=6      # number of separate fits to perform for each atom (more is better but needs more cores)
NTRY=3      # number of tries for each fit (more is better but slower)
NATOM=3     # total atoms in molecule to fit (usually number of atoms in molecule)
MAXCHG=4    # max charges per atom to fit (4 is usually enough for l=2 quality)

cd $WORKDIR
rm -r fit*

for ((j=1; j<=$NFIT; j++)); do
  cd $WORKDIR
  mkdir fit$j && cd fit$j

  for ((i=1; i<=$NATOM; i++)); do

    # Adapt this submission script for your cluster environment
    SCRIPT="atom"$i".sh"
    echo "#!/bin/bash

#SBATCH --job-name=atom$i
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=vshort

WORKDIR=$WORKDIR
BINDIR=$BINDIR
REFDIR=$REFDIR
MTPFILE=$MTPFILE
PCUBE=$PCUBE
DCUBE=$DCUBE
NTRY=$NTRY
ATOMINDEX=$i
MAXCHG=$MAXCHG

cd \$WORKDIR/fit$j
hostname

\$BINDIR/pcubefit.x -greedy -mtpfile \$MTPFILE -esp \$PCUBE -dens \$DCUBE -nacmin 1 -nacmax \$MAXCHG -atom \$ATOMINDEX -ntry \$NTRY -onlymultipoles -v > atom$i-fit.out
" > $SCRIPT

    sbatch $SCRIPT
  done
done

