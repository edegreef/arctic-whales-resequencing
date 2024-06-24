#!/bin/bash
#SBATCH --mail-user=evelien.degreef@umanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --time=12:00:00
#SBATCH --mem=80000M
#SBATCH --array=1-100
#SBATCH --job-name=estimate_timepoints_bow
#SBATCH --account=def-coling

# 100 array runs/iterations of 'smc++ estimate' for a population (this script is for one population)
# also remember to make the folders first ("smc_out_timepoints")

source /home/edegreef/smcpp/bin/activate

mkdir /scratch/edegreef/bowhead/smcpp_rightwhale/smc_out_timepoints/run_${SLURM_ARRAY_TASK_ID}/

echo smc_out_timepoints/run_${SLURM_ARRAY_TASK_ID}/

# Run smc++ estimate
# mutation rate 2.69E-08 
smc++ estimate 2.69e-8 -o /scratch/edegreef/bowhead/smcpp_rightwhale/smc_out_timepoints/run_${SLURM_ARRAY_TASK_ID}/ --regularization-penalty 4.0 --nonseg-cutoff 100000 --thinning 2000 --cores 4 --timepoints 10 10000 /scratch/edegreef/bowhead/smcpp_rightwhale/masked_data/ALL*

# Run again for other sample and other mu (1.20E-08)
#ALL.88_Pang
#ALL.WBF_2005_0298
