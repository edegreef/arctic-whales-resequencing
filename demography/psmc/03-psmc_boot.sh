#!/bin/bash

sample=BM_RB_2005_001_S43.rightwhale

cd BM_RB_2005_001/boot

# split psmcfa file for bootstrapping:
/home/degreefe/programs/psmc/utils/splitfa $sample.psmcfa > $sample.split.psmcfa

# run bootstrap x100 (i think -P is threads)
seq 100 | xargs -P 10 -i echo /home/degreefe/programs/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-{}.$sample.psmc $sample.split.psmcfa | sh

# combine with main run
cat $sample.N25.t15.r5.psmc round-*.$sample.psmc > combined.$sample.psmc

# see if plotting works (and then use the all the plot.#.txt files to plot in R, there would be 100 files per sample)
/home/degreefe/programs/psmc/utils/psmc_plot.pl -R -pY50000 combined combined.$sample.psmc
