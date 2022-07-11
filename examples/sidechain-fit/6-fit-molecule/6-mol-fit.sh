#!/bin/bash

# In this step we use the atomic fits from the previous step to create an initial
# guess, then we use differential evolution to fit a molecular charge model.
# Again, we repeat the fitting NTRY times for each run to get a good solution,
# and run NFIT separate fits to explore varied solutions.
#
# We also run fits with different numbers of charges, so that we can see the
# impact of adding more charges on RMSE and choose a good compromise between
# accuracy and computational cost
#
# The additional conformers are introduced at this stage, so that a model is
# found that optimally describes all provided conformers simultaneously.
#
# The file defining local reference axes is required to transform each 
# candidate charge model for each conformer used for fitting. The file contains
# the residue name, followed by groups of 3 bonded atoms with associated axis
# system preferred by the user (z-axis along a bond=BO, z-axis along the A-B-C
# bisector=BI) E.g.
# TRP
# 2 1 3 BO
# 4 3 5 BI
# ...
# Note that the choice of local axis type is not critical, all atoms should
# appear in at least one axis system, atom numbering corresponds to the
# reference CUBE files, and it's fine for an atom to appear in more than one
# local axis frame.

# Note: adapt the Slurm script below to whatever cluster environment you're using

# FOLDERS
ROOTDIR=/home/devereux/MDCM-release/examples/sidechain-fit
WORKDIR=$ROOTDIR/6-fit-molecule
BINDIR=$ROOTDIR/../../bin
REFDIR=$ROOTDIR/../ref
ATOMDIR=$ROOTDIR/5-fit-atoms
CUBEDIR=$ROOTDIR/3-modify-cubes
# INPUT
MTPFILE=$ROOTDIR/4-mtp-fit/fitted-mtpl.dat # from step 1
PCUBE=$CUBEDIR/equi.pot.mod.cube      # Equilibrium geometry ESP (frozen atoms removed)
DCUBE=$CUBEDIR/equi.dens.mod.cube     # Equilibrium geometry density (frozen atoms removed)
PCUBE2=$CUBEDIR/conf2.pot.mod.cube    # 2nd conformer ESP (frozen atoms removed)
DCUBE2=$CUBEDIR/conf2.dens.mod.cube   # 2nd conformer density (frozen atoms removed)
FRAMEFILE=$WORKDIR/trp-frames.txt     # file defining local reference axes
# FITTING PARAMETERS
NFIT=5        # number of separate fits to perform for each fragment (should not exceed atom fits)
NTRY=1        # number of tries for each fit
MAXATMCHG=4   # max charges per atom used in atom fits in step 2
MINCHG=18     # smallest number of charges to fit for molecule
MAXCHG=27     # largest number of charges to fit for molecule

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
PCUBE2=$PCUBE2
DCUBE2=$DCUBE2
FRAMEFILE=$FRAMEFILE
MTPFILE=$MTPFILE
NTRY=$NTRY
MINCHG=$NCHG
MAXCHG=$NCHG
NAME="${k}chgs"

cd $WORKDIR/fit$j

\$BINDIR/pcubefit.x -greedy -esp \$PCUBE -dens \$DCUBE -esp \$PCUBE2 -dens \$DCUBE2 -mtpfile \$MTPFILE -ncmin \$MINCHG -ncmax \$MAXCHG -nacmax 4 -ntry \$NTRY -frames \$FRAMEFILE -converge 0.01 -v > \$NAME".out"
" > $SCRIPT

    sbatch $SCRIPT
    sleep 1
  done
done

