#!/bin/bash

# In this step we divide the molecule into fragments and fit each fragment 
# separately If we create more fragments then we can better parallelize the
# fit, and each fit will be faster as fewer charges are required, but the error
# and redundancy will be larger when we combine them in the final molecular
# refinement step. With 18 charges per fragment expect several hours if default
# fitting parameters are used...
#
# Here we fit to 2 conformers simultaneously, reducing the conformational 
# sensitivity of the resulting models in subsequent MD simulations. This
# approach requires local reference axes to be defined using groups of 3 bonded
# atoms (not colinear) and supplied by the user with format:
#
# TRP
# 2 1 3 BO
# 4 3 5 BI
# ...
#
# The first line defines the residue name, the subsequent lines define the 
# atoms that comprise each local axis frame. Each atom must appear in at least
# one frame, if an atom appears in several frames then charges will be added in
# only the first. The label "BO" or "BI" requests that the z-axis lie long 
# either a bond or the A-B-C bisector, this option generally does not strongly
# affect fitting.

# Note: adapt the Slurm script below to whatever cluster environment you're using

# FOLDERS
ROOTDIR=/home/devereux/MDCM-release/examples/fragment-multiconf
WORKDIR=$ROOTDIR/3-fit-frags
BINDIR=$ROOTDIR/../../bin
REFDIR=$ROOTDIR/ref
ATOMDIR=$ROOTDIR/2-fit-atoms
# INPUT
MTPFILE=$ROOTDIR/1-mtp-fit/fitted-mtpl-equi.dat # from step 1
PCUBE=$REFDIR/equi.pot.mod.cube    # not really used as we fit to fragment ESP
DCUBE=$REFDIR/equi.dens.mod.cube   
MTPFILE2=$ROOTDIR/1-mtp-fit/fitted-mtpl-conf2.dat # from step 1
PCUBE2=$REFDIR/conf2.pot.mod.cube  # not really used as we fit to fragment ESP
DCUBE2=$REFDIR/conf2.dens.mod.cube
FRAMEFILE=$REFDIR/trp-frames.txt   # this reference axis file is supplied by the user
# user-defined relaxed cutoff tolerance for early exit from DE (optional, faster):
CONV="0.01" 

# FITTING PARAMETERS
NFIT=5        # number of separate fits to perform for each fragment (should not exceed atom fits)
NTRY=2        # number of tries for each fit
MAXATMCHG=4   # max charges per atom used in atom fits
NFRAG=2       # number of fragments to fit
MINFRAGCHG=8  # smallest number of charges to fit for fragment
MAXFRAGCHG=20 # largest number of charges to fit for fragment
ATOMLIST=("1,2,3,4,5,6,7,8" "9,10,11,12,13,14,15,16,17,18") # array of atom indices of each fragment

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
#SBATCH --ntasks=8
#SBATCH --partition=infinite

WORKDIR=$WORKDIR
BINDIR=$BINDIR
REFDIR=$REFDIR
PCUBE=$PCUBE
DCUBE=$DCUBE
MTPFILE=$MTPFILE
PCUBE2=$PCUBE2
DCUBE2=$DCUBE2
MTPFILE2=$MTPFILE2
NTRY=$NTRY
ATOMLIST=\"${ATOMLIST[$ii]}\"
MINCHG=$NCHG
MAXCHG=$NCHG
NAME="${k}chgs"
FRAMEFILE=$FRAMEFILE
CONV="$CONV"

cd $WORKDIR/frag$i/fit$j

\$BINDIR/pcubefit.x -greedy -esp \$PCUBE -dens \$DCUBE -mtpfile \$MTPFILE -esp \$PCUBE2 -dens \$DCUBE2 -mtpfile \$MTPFILE2 -ncmin \$MINCHG -ncmax \$MAXCHG -atom \$ATOMLIST -nacmax 4 -ntry \$NTRY -frames \$FRAMEFILE -converge \$CONV -v > \$NAME".out"
" > $SCRIPT

      sbatch $SCRIPT
      sleep 1
    done
  done
done


