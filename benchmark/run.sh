#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=medium
#SBATCH --mem=0

ml lang JuliaHPC

respath="results"
rm -rf $respath
threads=(1 4 8)
for t in ${threads[@]}; do
    julia --project=.. -t $t bbvv.jl $respath
done
