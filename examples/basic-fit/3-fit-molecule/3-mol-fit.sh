#!/bin/bash

# In this step we use the atomic fits from the previous step to create an initial
# guess, then we use differential evolution to fit a molecular charge model.
# Again, we repeat the fitting NTRY times for each run to get a good solution,
# and run NFIT separate fits to explore varied solutions.

# We can also run fits with different numbers of charges, so that we can see the
# impact of adding more charges on RMSE and choose a good compromise between
# accuracy and computational cost

# Note: adapt the Slurm script below to whatever cluster environment you're using

# FOLDERS
ROOTDIR=/home/devereux/MDCM-release/examples/basic-fit
WORKDIR=$ROOTDIR/3-fit-molecule
BINDIR=$ROOTDIR/../../bin
REFDIR=$ROOTDIR/../ref
ATOMDIR=$ROOTDIR/2-fit-atoms
# INPUT
MTPFILE=$ROOTDIR/1-mtp-fit/fitted-mtpl.dat # from step 1
PCUBE=$REFDIR/h2o-pot.cube
DCUBE=$REFDIR/h2o-dens.cube

# FITTING PARAMETERS
NFIT=6        # number of separate fits to perform for each fragment (should not exceed atom fits)
NTRY=3        # number of tries for each fit
MAXATMCHG=4   # max charges per atom used in atom fits in step 2
MINCHG=4      # smallest number of charges to fit for molecule
MAXCHG=10     # largest number of charges to fit for molecule

cd $WORKDIR

for ((j=1; j<=$NFIT; j++)); do
  cd $WORKDIR
  mkdir -p fit$j && cd fit$j
  cp $ATOMDIR/fit$j/*xyz .  # copy atom fit results used to build starting population

  for ((k=$MINCHG; k<=$MAXCHG; k++)); do
    SCRIPT=$k"chgs.sh"
    NCHG=$k
    echo "#!/bin/bash

#SBATCH --job-name=fit${j}q${k}
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=long

WORKDIR=$WORKDIR
BINDIR=$BINDIR
REFDIR=$REFDIR
PCUBE=$PCUBE
DCUBE=$DCUBE
MTPFILE=$MTPFILE
NTRY=$NTRY
MINCHG=$NCHG
MAXCHG=$NCHG
NAME="${k}chgs"

cd $WORKDIR/fit$j

\$BINDIR/pcubefit.x -greedy -esp \$PCUBE -dens \$DCUBE -mtpfile \$MTPFILE -ncmin \$MINCHG -ncmax \$MAXCHG -nacmax 4 -ntry \$NTRY -v > \$NAME".out"
" > $SCRIPT

    sbatch $SCRIPT
    sleep 1
  done
done

