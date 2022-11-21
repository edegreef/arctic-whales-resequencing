#!/bin/bash

source /home/edegreef/smcpp/bin/activate

smc++ plot -c smc_1group_plot.pdf smc_out_timepoints/1group/run_*/model.final.json
