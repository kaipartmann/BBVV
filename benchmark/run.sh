#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=medium
#SBATCH --mem=0
#SBATCH --reservation=hackathon

ml lang JuliaHPC

respath="results_150"
rm -rf $respath
threads=(1 2 4 8 16 32 64)
for t in ${threads[@]}; do
    julia --project=.. -t $t bbvv.jl $respath
done
