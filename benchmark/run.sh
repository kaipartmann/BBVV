#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=kai.partmann@uni-siegen.de
# # SBATCH --reservation=hackathon
# # SBATCH --account=hpc-lco-usrtr
#SBATCH --exclusive # good for benchmarking!

ml lang JuliaHPC

respath="results_150"
rm -rf $respath
threads=(1 2 4 8 16 32 64)
for t in ${threads[@]}; do
    julia --project=.. -t $t bbvv.jl $respath
done
