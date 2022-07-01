#!/bin/bash

# In this step we run a symmetry-constrained fit for the whole molecule. 5 separate fits
# are performed here, for each fit we try fitting models with between MINMOLCHG and
# MAXMOLCHG charges in total. For some numbers of charges no model is possible
# without breaking symmetry constraints, so these models are skipped.

# Note: adapt the Slurm script below to whatever cluster environment you're using

# FOLDERS
ROOTDIR=/pchem-data/meuwly/devereux/MDCM/examples/benzene-sym
WORKDIR=$ROOTDIR/4-fit-molecule
BINDIR=/pchem-data/meuwly/devereux/MDCM/bin
REFDIR=$ROOTDIR/ref
ATOMDIR=$ROOTDIR/3-fit-atoms
# INPUT
MTPFILE=$ROOTDIR/2-mtp-fit/fitted-mtpl-l5.dat # from step 2
PCUBE=$REFDIR/benzene-pot.cube # again, not really used for fragments
DCUBE=$REFDIR/benzene-dens.cube
# FITTING PARAMETERS
NFIT=5       # number of separate fits to perform for each fragment (should not exceed atom fits)
NTRY=4        # number of tries for each fit
MAXATMCHG=4   # max charges per atom used in atom fits
MINMOLCHG=12   # smallest number of charges to fit for fragment
MAXMOLCHG=36 # largest number of charges to fit for fragment

cd $WORKDIR

for ((j=1; j<=$NFIT; j++)); do
  cd $WORKDIR
  mkdir -p fit$j && cd fit$j
  cp $ATOMDIR/fit$j/*xyz .  # copy atom fit results used to build starting population
  cp $ATOMDIR/fit$j/*fit .  # and fitted symmetry parameters

  for ((k=$MINMOLCHG; k<=$MAXMOLCHG; k++)); do
    SCRIPT=$k"chgs.sh"
    NCHG=$k
    echo "#!/bin/bash

#SBATCH --job-name=frag${i}fit${j}q${k}
#SBATCH --nodes=1
#SBATCH --ntasks=4
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
MINCHG=$NCHG
MAXCHG=$NCHG
NAME="${k}chgs"

cd $WORKDIR/fit$j
echo frag${i}fit${j}q${k}
hostname

module load gcc/gcc-7.4.0-openmpi-3.1.4-cuda-for-11.4

\$BINDIR/cudacubefit.x -greedy \$MTPFILE -esp \$PCUBE -dens \$DCUBE -ncmin \$MINCHG -ncmax \$MAXCHG -nacmax 4 -ntry \$NTRY -sym -v -gpu > \$NAME".out"
" > $SCRIPT

    sbatch $SCRIPT
    sleep 1
  done
done

