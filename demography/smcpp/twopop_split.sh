#!/bin/bash

#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --time=2:00:00
#SBATCH --mem=80000M
#SBATCH --array=1-20
#SBATCH --job-name=twopop_split

source ~/smcpp/bin/activate

pops=ABGF
pop1=AB
pop2=GF
masked_data_folder=masked_data_run_ABGF

mkdir split_estimate
mkdir split_estimate/ABGF

echo run_${SLURM_ARRAY_TASK_ID}

smc++ split -o split_estimate/$pops/combined_run${SLURM_ARRAY_TASK_ID} smc_out_timepoints/$pop1/run_${SLURM_ARRAY_TASK_ID}/model.final.json smc_out_timepoints/$pop2/run_${SLURM_ARRAY_TASK_ID}/model.final.json $masked_data_folder/*.smc
