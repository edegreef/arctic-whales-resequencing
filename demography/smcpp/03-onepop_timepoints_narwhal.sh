#!/bin/bash
#SBATCH --mail-user=evelien.degreef@umanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --time=18:00:00
#SBATCH --mem=80000M
#SBATCH --array=1-100
#SBATCH --job-name=estimate_timepoints_nar

# 100 array runs/iterations of 'smc++ estimate' for a population (this script is for one population)
# also remember to make the folders first ("smc_out_timepoints", and a folder for each pop ID if running more later (here I used "1group"))

source /home/edegreef/smcpp/bin/activate

mkdir /scratch/edegreef/whales/NAR_smcpp/smc_out_timepoints/1group_2D2/run_${SLURM_ARRAY_TASK_ID}/

echo smc_out_timepoints/1group_2D2/run_${SLURM_ARRAY_TASK_ID}/

# Run smc++ estimate
smc++ estimate 1.56e-8 -o /scratch/edegreef/whales/NAR_smcpp/smc_out_timepoints/1group_2D2/run_${SLURM_ARRAY_TASK_ID}/ --regularization-penalty 4.0 --nonseg-cutoff 100000 --thinning 2000 --cores 4 --timepoints 10 10000 /scratch/edegreef/whales/NAR_smcpp/masked_data/ALL.*
