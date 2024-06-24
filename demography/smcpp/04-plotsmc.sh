#!/bin/bash

# here just quick plots, but use the smcpp output files in R to make better plots

source /home/edegreef/smcpp/bin/activate

smc++ plot -c smc_bowhead_RWmap_2.69e-8_plot.pdf smc_out_timepoints/run_*/model.final.json
smc++ plot -c smc_narwhal_EBB_1.56e-8_plot.pdf smc_out_timepoints_subgroup/EBB/run_*/model.final.json
smc++ plot -c smc_narwhal_RB_1.56e-8_plot.pdf smc_out_timepoints_subgroup/RB/run_*/model.final.json
smc++ plot -c smc_narwhal_WBB_1.56e-8_plot.pdf smc_out_timepoints_subgroup/WBB/run_*/model.final.json

# other mus
smc++ plot -c smc_bowhead_RWmap_1.20e-8_plot.pdf smc_out_timepoints_othermu/run_*/model.final.json
smc++ plot -c smc_narwhal_EBB_1.41e-8_plot.pdf smc_out_timepoints_subgroup_othermu/EBB/run_*/model.final.json
smc++ plot -c smc_narwhal_RB_1.41e-8_plot.pdf smc_out_timepoints_subgroup_othermu/RB/run_*/model.final.json
smc++ plot -c smc_narwhal_WBB_1.41e-8_plot.pdf smc_out_timepoints_subgroup_othermu/WBB/run_*/model.final.json
