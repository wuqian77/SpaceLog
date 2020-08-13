#!/usr/bin/env
#!/bin/bash
#$ -cwd

#-- Important directories --------------------------------------
pg_dir=/home/wu_v/
lg_dir=/home/wu_v/log

file=step_1_GTex_gene_expression_byGeneSet_short.R
cd $pg_dir

#-- Submit jobs --------------------------------------
for k in `seq 1 189`
do
sbatch --wrap="R CMD BATCH --no-save --no-restore '--args $k ' $file $lg_dir/$file.$k.out"
done





