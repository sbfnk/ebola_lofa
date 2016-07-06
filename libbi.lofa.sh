#!/bin/bash
#$ -N lofa
#$ -V -cwd
#$ -M sebastian.funk@lshtm.ac.uk
#$ -l mem_free=4000M,h_vmem=8000M
#$ -pe smp 8
#$ -q parallel.q
#$ -R y
#$ -t 1-20

r0_trajectory=$1
if [ -z "$r0_trajectory" ]
then
    echo Error: no r0 trajectory given
    exit 1
fi

seed=$RANDOM

echo Seed: $seed

source ~/.bashrc
source ~/perl5/perlbrew/etc/bashrc

##Rscript ~/code/ebola_lofa/R/sample_posterior.r -n 10000 -p 1000 -r $r0_trajectory -t 4 -e $seed -i 10 -q -l -f -k -o ebola_lofa_poisson
Rscript ~/code/ebola_lofa/R/sample_posterior.r -n 10000 -p 1000 -r $r0_trajectory -t 8 -e $seed -i 10 -q -l -f -o ebola_lofa_$r0_trajectory\_poisson_${SGE_TASK_ID} -c 16384 -b ${SGE_TASK_ID} 
