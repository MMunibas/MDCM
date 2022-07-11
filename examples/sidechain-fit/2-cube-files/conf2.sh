#!/bin/bash

#SBATCH --job-name=conf2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --cpus-per-task=8
# 
hostname

 
#set -xv
export g16root=/opt/cluster/programs/g16-c.01
export GAUSS_SCRDIR=/scratch/devereux/conf2
mkdir -p /scratch/devereux/conf2
source $g16root/g16/bsd/g16.profile

$g16root/g16/g16 /home/devereux/MDCM-release/examples/sidechain-fit/2-cube-files/conf2.com /scratch/devereux/conf2/conf2.out

# don't delete the result file if not able to copy to fileserver 
cp /scratch/devereux/conf2/conf2.out /home/devereux/MDCM-release/examples/sidechain-fit/2-cube-files/conf2.out 
status=$?
if [ $status -eq 0 ] 
then 
   rm -rf /scratch/devereux/conf2
else
   host=`/bin/hostname`
   /usr/bin/Mail -v -s "Error at end of batch job" $USER <<EOF

At the end of the batch job the system could not copy the output file
	$host:/scratch/devereux/conf2/conf2.out
to
	/home/devereux/MDCM-release/examples/sidechain-fit/2-cube-files/conf2.out
Please copy this file by hand or inform the system manager.

EOF
 
fi
