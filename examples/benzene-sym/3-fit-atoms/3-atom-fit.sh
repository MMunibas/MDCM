#!/bin/bash

# Here we fit symmetry-constrained atomic MDCM models to the atomic multipoles from
# the previous step using a GPU. Only one atom of each symmetry-related group is
# fitted. To improve quality we should run more fits to choose the best, and also
# increase the number of DE tries for each fit. Here we show an example with 5 fits
# and 8 refinement attempts (NTRY) per fit. Note that multiple fits are more
# important when running with symmetry-constraints as different constraints may be
# automatically selected for each fit.

# Note: adapt the Slurm script below to whatever cluster you're using...

# FOLDERS
ROOT=/pchem-data/meuwly/devereux/MDCM/examples/benzene-sym
WORKDIR=$ROOT/3-fit-atoms
BINDIR=/pchem-data/meuwly/devereux/MDCM/bin
REFDIR=$ROOT/ref
MTPDIR=$ROOT/2-mtp-fit
# INPUT
PCUBE=$REFDIR/benzene-pot.cube # not really used here as we fit to multipolar ESP
DCUBE=$REFDIR/benzene-dens.cube
MTPFILE=$MTPDIR/fitted-mtpl-l5.dat # from step 2
# FITTING PARAMETERS
NFIT=5    # number of separate fits to perform for each atom (more is better but needs more cores)
NTRY=8      # number of tries for each fit (more is better but slower)
NATOM=12     # total atoms in molecule to fit (usually number of atoms in molecule)
MAXCHG=4      # max charges per atom to fit (4 is usually enough for l=2 quality)

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
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint=RTX2080Ti

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

module load gcc/gcc-7.4.0-openmpi-3.1.4-cuda-for-11.4

\$BINDIR/cudacubefit.x -greedy \$MTPFILE -esp \$PCUBE -dens \$DCUBE -nacmin 1 -nacmax \$MAXCHG -atom \$ATOMINDEX -ntry \$NTRY -onlymultipoles -sym -gpu -v > atom$i-fit.out
" > $SCRIPT

    sbatch $SCRIPT
  done
done

