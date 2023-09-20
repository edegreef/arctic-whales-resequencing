#!/bin/bash

# run GONE in the GONE_program folder (downloaded files from https://github.com/esrud/GONE)

cd GONE_program

# run the script_GONE.sh with parameter file (INPUT_PARAMETERS_FILE)
./script_GONE.sh narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n57.scafrename

# run bowhead one after narwhal is done & moved to separate folder
