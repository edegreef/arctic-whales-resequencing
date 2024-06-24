#!/bin/bash

# convert vcf to format for GONE
# ./GONE_prep.sh narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.EBB

# using vcfs that are already downsampled to population-level

# script:
snps=$1
scaf_list=narwhal_10MB_rename_scaffolds.txt
chr=21 #21 for narwhal, 20 for bowhead

# Need to name each scaffold from 1-# to make format compatible with GONE.
# run these first
#grep -v "^#" bowhead_RWmap_snps_rmfix.filter1.miss.biallel.min100kb.autosomes.hwe.n20.10MB.vcf | cut -f1 | uniq > bowhead_RWmap_scaf_10MB_list_check.txt
#grep -v "^#" narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.EBB.10MB.vcf | cut -f1 | uniq > narwhal_scaf_10MB_list_check.txt

# add row number to become the new scaffold names:
#awk '{print $0,NR}' bowhead_RWmap_scaf_10MB_list_check.txt > bowhead_RWmap_10MB_rename_scaffolds.txt
#awk '{print $0,NR}' narwhal_scaf_10MB_list_check.txt > narwhal_10MB_rename_scaffolds.txt


# Rename scaffolds-- in hindsight, better do to this step before separating vcfs into each pop.
/home/degreefe/programs/bcftools-1.9/bcftools annotate --rename-chrs $scaf_list $1.vcf -Oz -o $1.scafrename.vcf.gz

# Convert to ped/map format (*important to do this with plink and not vcftools)
/home/degreefe/programs/plink --vcf $1.scafrename.vcf.gz --recode --out $1.scafrename --double-id --chr-set $chr
