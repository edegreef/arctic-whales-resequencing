#!/bin/bash
#SBATCH --mail-user=evelien.degreef@umanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --time=12:00:00
#SBATCH --mem=80000M
#SBATCH --array=1-100
#SBATCH --job-name=estimate_timepoints_narwhal_sub
#SBATCH --account=def-coling

# 100 array runs/iterations of 'smc++ estimate' for a population (this script is for one population)
# also remember to make the folders first ("smc_out_timepoints_subgroup", and a folder for each pop)

source /home/edegreef/smcpp/bin/activate

# Run smc++ estimate

# Canadiah High Arctic / Western Baffin Bay
mkdir /scratch/edegreef/narwhal_smcpp/smc_out_timepoints_subgroup/WBB/run_${SLURM_ARRAY_TASK_ID}/
echo smc_out_timepoints_subgroup/WBB/run_${SLURM_ARRAY_TASK_ID}/
# Run smc++ estimate
smc++ estimate 1.56e-8 -o /scratch/edegreef/narwhal_smcpp/smc_out_timepoints_subgroup/WBB/run_${SLURM_ARRAY_TASK_ID}/ --regularization-penalty 4.0 --nonseg-cutoff 100000 --thinning 2000 --cores 4 --timepoints 10 10000 /scratch/edegreef/narwhal_smcpp/masked_data_subgroup/WBB.*

# Baffin Island / Eastern Baffin Bay
mkdir /scratch/edegreef/narwhal_smcpp/smc_out_timepoints_subgroup/EBB/run_${SLURM_ARRAY_TASK_ID}/
echo smc_out_timepoints_subgroup/EBB/run_${SLURM_ARRAY_TASK_ID}/
# Run smc++ estimate
smc++ estimate 1.56e-8 -o /scratch/edegreef/narwhal_smcpp/smc_out_timepoints_subgroup/EBB/run_${SLURM_ARRAY_TASK_ID}/ --regularization-penalty 4.0 --nonseg-cutoff 100000 --thinning 2000 --cores 4 --timepoints 10 10000 /scratch/edegreef/narwhal_smcpp/masked_data_subgroup/EBB.*

# Repulse Bay / northern Hudson Bay
mkdir /scratch/edegreef/narwhal_smcpp/smc_out_timepoints_subgroup/RB/run_${SLURM_ARRAY_TASK_ID}/
echo smc_out_timepoints_subgroup/RB/run_${SLURM_ARRAY_TASK_ID}/
# Run smc++ estimate
smc++ estimate 1.56e-8 -o /scratch/edegreef/narwhal_smcpp/smc_out_timepoints_subgroup/RB/run_${SLURM_ARRAY_TASK_ID}/ --regularization-penalty 4.0 --nonseg-cutoff 100000 --thinning 2000 --cores 4 --timepoints 10 10000 /scratch/edegreef/narwhal_smcpp/masked_data_subgroup/RB.*
