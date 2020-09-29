#!/bin/bash

# Script to find the best combination of fragment results that combine to a given number
# of molecule charges

# TODO: slow for larger systems with many permutations - should rather read all
# RMSE values into memory at start and find best RMSE using e.g. Python rather than BASH

# FOLDERS
ROOT=/home/devereux/MDCM-git/examples/naptha
WORKDIR=$ROOT/5-combine-frags
BINDIR=/home/devereux/MDCM-git/bin
FRAGDIR=$ROOT/4-fit-frags # folder containing fitted fragments
# FITTING PARAMETERS
NFRAG=5             #number of fragments fitted
NFIT=1              #number of fits per fragment
MINCHGS=18          #minimum number of charges for whole molecule
MAXCHGS=64          #maximum number of charges for whole molecule
MINFRAGCHGS=6       #minimum charges used to fit fragments in previous step
MAXFRAGCHGS=14      #maximum charges used to fit fragments in previous step

cd $FRAGDIR

for ((i=$MINCHGS; i<=$MAXCHGS; i++)); do

  # get all permutations that add up to $i using $NFRAG terms
  $BINDIR/charge-permutations.py $NFRAG $i $MINFRAGCHGS $MAXFRAGCHGS > permutations-${i}.dat
  nperm=$(sed -n $= permutations-${i}.dat)
  [ -z $nperm ] && echo "Warning, no solution for $i charges" && continue
  best_rmse=999999
  for ((j=1; j<=$nperm; j++)); do
    perm=$(sed -n ${j}p permutations-${i}.dat)
    for ((k=1; k<=$NFRAG; k++)); do
      nchg[$k]=$(echo $perm | awk "{print \$${k}}")
      rmse[$k]=999999
      for ((l=1; l<=$NFIT; l++)); do
        if [ -e $FRAGDIR/frag$k/fit$l/${nchg[$k]}charges.xyz ]; then
          trmse=$(grep RMSE $FRAGDIR/frag$k/fit$l/${nchg[$k]}charges.xyz | awk '{print $2}' | sed "s/E/*10^/g" | sed "s/+//g")
          if (( $(echo "$trmse < ${rmse[$k]}" | bc -l) )); then
            rmse[$k]=$trmse
            best[$k]=$l
            fragfile[$k]=$FRAGDIR/frag$k/fit$l/${nchg[$k]}charges.xyz
          fi
        fi
      done
    done
    fraglist=""
    mean_rmse=0
    for ((k=1; k<=$NFRAG; k++)); do
      mean_rmse=$(echo "$mean_rmse + ${rmse[$k]}" | bc -l)
      fraglist="$fraglist ${fragfile[$k]}"
    done
    mean_rmse=$(echo "$mean_rmse / $NFRAG" | bc -l)
    echo "$i chgs permutation $perm rmse=$mean_rmse"
    if (( $(echo "$mean_rmse < $best_rmse" | bc -l) )); then
      best_rmse=$mean_rmse
      echo "NEW BEST $i charges: $perm rmse = $mean_rmse"
      $BINDIR/combine-frags.py $fraglist
      echo >> combined.xyz
      echo "total molecule built from fragment files: $fraglist" >> combined.xyz
      mv combined.xyz $WORKDIR/$i-combined.xyz
    fi
  done
  rm permutations-${i}.dat

done
