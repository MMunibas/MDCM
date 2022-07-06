#!/bin/bash

# Here we fit atomic MDCM models to the atomic multipoles from the previous 
# step. To improve quality we should run more fits to choose the best, and also
# increase the number of DE tries for each fit. Here we show an example with 5 
# fits and 4 refinement attempts (NTRY) per fit.

# Note: adapt the Slurm script below to whatever cluster you're using...

# FOLDERS
ROOT=/home/devereux/MDCM-release/examples/fragment-fit
WORKDIR=$ROOT/2-fit-atoms
BINDIR=/home/devereux/MDCM-release/bin
REFDIR=$ROOT/../ref
MTPDIR=$ROOT/1-mtp-fit
# INPUT
PCUBE=$REFDIR/benzene-pot.cube # not really used here as we fit to multipolar ESP
DCUBE=$REFDIR/benzene-dens.cube
MTPFILE=$MTPDIR/fitted-mtpl.dat # from step 2
# FITTING PARAMETERS
NFIT=5    # number of separate fits to perform for each atom (more is better but needs more cores)
NTRY=4      # number of tries for each fit (more is better but slower)
NATOM=12     # total atoms in molecule to fit (usually number of atoms in molecule)
MAXCHG=4      # max charges per atom to fit (4 is usually enough for l=2 multipolar quality)

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
#SBATCH --ntasks=4
#SBATCH --partition=long

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

