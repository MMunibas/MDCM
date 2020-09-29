#!/bin/bash

# In this step we divide the molecule into fragments and fit each fragment separately
# If we create more fragments then we can better parallelize the fit, and each fit
# will be faster as fewer charges are required, but the error and redundancy will be
# larger when we combine them in the final molecular refinement step. With 18 charges
# per fragment expect several hours if default fitting parameters are used...

# Note: adapt the Slurm script below to whatever cluster environment you're using

# FOLDERS
ROOTDIR=/home/devereux/MDCM-git/examples/naptha/
WORKDIR=$ROOTDIR/4-fit-frags
BINDIR=/home/devereux/MDCM-git/bin
REFDIR=$ROOTDIR/ref
ATOMDIR=$ROOTDIR/3-fit-atoms
# INPUT
MTPFILE=$ROOTDIR/2-mtp-fit/fitted-mtpl.dat # from step 2
PCUBE=$ROOTDIR/2-mtp-fit/ditriantapole_expansion.cube # again, not really used
DCUBE=$REFDIR/naphta04.dens.cube
# FITTING PARAMETERS
NFIT=1        # number of separate fits to perform for each fragment (should not exceed atom fits)
NTRY=3        # number of tries for each fit
MAXATMCHG=4   # max charges per atom used in atom fits
NFRAG=5       # number of fragments to fit
MINFRAGCHG=6  # smallest number of charges to fit for fragment
MAXFRAGCHG=14 # largest number of charges to fit for fragment
ATOMLIST=("1,2,11,12" "5,8,14,16" "4,7,13,15" "9,10,17,18" "3,6") # array of atom indices of each fragment

cd $WORKDIR

for ((i=1; i<=$NFRAG; i++)); do
  cd $WORKDIR
  ii=$(($i-1))
  mkdir -p "frag"$i && cd "frag"$i

  for ((j=1; j<=$NFIT; j++)); do
    cd $WORKDIR"/frag"$i
    mkdir -p fit$j && cd fit$j
    cp $ATOMDIR/fit$j/*xyz .  # copy atom fit results used to build starting population

    for ((k=$MINFRAGCHG; k<=$MAXFRAGCHG; k++)); do
      SCRIPT=$k"chgs.sh"
      NCHG=$k
      echo "#!/bin/bash

#SBATCH --job-name=frag${i}fit${j}q${k}
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
ATOMLIST=\"${ATOMLIST[$ii]}\"
MINCHG=$NCHG
MAXCHG=$NCHG
NAME="${k}chgs"

cd $WORKDIR/frag$i/fit$j

\$BINDIR/pcubefit.x -greedy \$MTPFILE -esp \$PCUBE -dens \$DCUBE -ncmin \$MINCHG -ncmax \$MAXCHG -atom \$ATOMLIST -nacmax 4 -ntry \$NTRY -v > \$NAME".out"
" > $SCRIPT

      sbatch $SCRIPT
      sleep 1
    done
  done
done

