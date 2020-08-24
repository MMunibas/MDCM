#!/bin/bash

# Script to find the best combination of fragment results that combine to a given number
# of molecule charges

WORKDIR=/home/devereux/naphta04/mdcm_fit
BINDIR=/home/devereux/mdcm-fitting-code/mdcm-code
FRAGDIR=$WORKDIR/frag-fit
NFRAG=3             #number of fragments fitted
NFIT=2             #number of fits per fragment
MINCHGS=25          #minimum number of charges for whole molecule
MAXCHGS=36          #maximum number of charges for whole molecule

cd $FRAGDIR

for ((i=$MINCHGS; i<=$MAXCHGS; i++)); do

  # get all permutations that add up to $i using $NFRAG terms
  $BINDIR/charge-permutations.py $NFRAG $i > permutations-${i}.dat
  nperm=$(sed -n $= permutations-${i}.dat)
  best_rmse=999999
  for ((j=1; j<=$nperm; j++)); do
    perm=$(sed -n ${j}p permutations-${i}.dat)
    exists=true
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
      if [ "${rmse[$k]}" == "999999" ]; then
        exists=false  # no fit was performed with this many charges for this fragment
      fi
    done
    if $exists; then
      fraglist=""
      mean_rmse=0
      for ((k=1; k<=$NFRAG; k++)); do
        mean_rmse=$(echo "$mean_rmse + ${rmse[$k]}" | bc -l)
        fraglist="$fraglist ${fragfile[$k]}"
      done
      mean_rmse=$(echo "$mean_rmse / $NFRAG" | bc -l)
      if (( $(echo "$mean_rmse < $best_rmse" | bc -l) )); then
        best_rmse=$mean_rmse
        echo "NEW BEST $i charges: $perm rmse = $mean_rmse"
        $BINDIR/combine-frags.py $fraglist
        echo >> combined.xyz
        echo "total molecule built from fragment files: $fraglist" >> combined.xyz
        mv combined.xyz $i-combined.xyz
      fi
    fi
  done
  rm permutations-${i}.dat

done
