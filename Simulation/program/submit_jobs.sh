#!/bin/bash
#$ -cwd

#-- Important directories --------------------------------------
pg_dir=Program
lg_dir=Program/log
#-- Submit jobs --------------------------------------
file=step_1_1_SpaceLog.R
step=simu
queue=A
job=simu
cd $pg_dir
n=400

for method in `seq 2 4`
do
for k in `seq 1 100`
do
for p in 100 200 300 
do
for e in `seq 1 2`
do
for model in `seq 1 2`
do
sbatch --wrap="R CMD BATCH --no-save --no-restore '--args $method $k $n $p $e $model' $file $lg_dir/$file.$method.$k.$n.$p.$e.$model.out"
done
done
done
done
done


