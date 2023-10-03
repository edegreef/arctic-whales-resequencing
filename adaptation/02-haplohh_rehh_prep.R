#!/usr/bin/Rscript

# Run data2haplohh in rehh to create the wgscan file
# Rehh vignette: https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html
# Help from Matt Thorstensen to set this up


# Ran this step on the biology cluster because going through all scaffolds and takes time

# Load modules
library(rehh)
library(vcfR)
library(tidyverse)

###############################################################
###### Step 1: Load vcf and create map file
###############################################################

setwd("/home/degreefe/whales/narwhal/rehh")

vcf <- read.vcfR("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.n57.mac2.miss01.imputed.vcf.gz")
chroms <- unique(vcf@fix[,1])

# make map file, has 5 columns
map_file <- tibble(marker_ID = paste(vcf@fix[,1],"_",vcf@fix[,2],sep=""),
                   chromosome = vcf@fix[,1], 
                   position = vcf@fix[,2], 
                   ref = vcf@fix[,4], 
                   alt = vcf@fix[,5])

# Save map file into folder
write_delim(map_file, "mapfile_narwhal.inp", delim = "\t", col_names = FALSE)



#### Step 2 actually run the haplohh

for(i in 1:length(chroms)) {
  
  # haplotype file name for each scaffold
  hap_file = "narwhal_snps.filter1.miss.biallel.min100kb.autosomes.n57.mac2.miss01.imputed.vcf.gz"
  
  #create internal representation
  hh <- data2haplohh(hap_file = hap_file,
                     map_file = "mapfile_narwhal.inp", 
                     chr.name = chroms[i],
                     polarize_vcf = FALSE,
                     vcf_reader = "vcfR")
  
  # perform scan on a single chromosome (calculate iHH values)
  scan <- scan_hh(hh, polarized = FALSE)
  
  # concatenate chromosome-wise data frames to a data frame for the whole genome
  if (i == 1) {
    wgscan <- scan
  } else {
    wgscan <- rbind(wgscan, scan)
  }
}

# save as file
write.csv(wgscan, "narwhal_allchr_wgscan.csv")

# narwhal starting on Jan 29 3:15pm

# save object
saveRDS(wgscan, file="narwhal_wgscan.rds")


# Repeat for bowhead whale file
