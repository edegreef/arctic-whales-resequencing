# Script for preparing the adaptive SNP dataset to use in Gradient Forest

# Main steps:
### 1) Run genome scan with SNMF, then 2) run LFMMs with environmental predictors, and then 3) prepare snp list to extract to use for Gradient Forest

# Make sure to use snps that have been filtered for MAF, and thinned/pruned. starting with ped/map format for snmf

# Load packages
library(LEA)
library(vcfR)
library(tidyverse)
library(qvalue)

# LEA doesn't work well in a dropbox folder so doing on desktop for now
setwd("C:/Users/eveli/Desktop/LEA_narwhal_GF_rerun")

########## PART 1: SNMF
# Convert ped file to geno file (automatically outputs in working directory)
ped2geno("narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.ped")

# Run snmf for K 1-4, if want to save time can just use a specific K based on previous snmf run and/or fewer repetitions.
project=NULL
project=snmf("narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.geno", K=1:4, entropy=TRUE, repetitions=10, project="new")

# Can also load snmf project if saved
project=load.snmfProject("narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.snmfProject")

# Genome scan for selection: population differentiation tests. can't do this on K=1. Use the best K for the dataset. Here using K=2 for narwhal.
p <- snmf.pvalues(project, entropy=TRUE, ploidy=2, K=2)
pvalues <- p$pvalues

# plot histogram of pvalues
par(mfrow=c(1,1))
hist(pvalues, col="orange")
plot(-log10(pvalues), pch=19, col="blue", cex=.5)

# just to look at data frame and log vals
pvals <- as.data.frame(pvalues)
pvals$log <- -log10(pvals$pvalues)

# how many with -log10 minimum 3?
pvals_log10_3 <- subset(pvals, log >= 3)

# Because not many significant snps in the narwhal genome scan, going to use top 1%

# Adding snp info to the p values
# load map file
snp_info <- read.table("narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.map")

# pull out only CHROM and POS
snp_info <- select(snp_info, V1, V4)

# get the log p's, add snp_info, order from large to small, pull top 10%
logp <- as.data.frame(-log10(pvalues))
logp$original_order <- 1:nrow(logp)
snps_logp <- cbind(snp_info, logp)
colnames(snps_logp) <- c("CHROM", "POS", "logp", "order")

snps_logp2 <- snps_logp[,c(4,1,2,3)]

# putting in order of logp (higher more significant)
snps_logp_order <- snps_logp2[order(-snps_logp2$logp),]

# set variable for number of snpts
count <- nrow(snps_logp_order)

# extract top x%, using 1% here
top_count <- floor(0.01*count)
snp_subset <- as.data.frame(snps_logp_order[1:top_count,])

# add line (smallest log value 1.86 when looking at snp_subset) to plot just to see
abline(h=1.86, col="red")


# save list of snps - base genome scan, save top snps as .txt or .csv
write_delim(snp_subset, "narwhal_K2.imputed.thin1000.top0.01logp.subset.txt")
write.csv(snp_subset, "narwhal_K2.imputed.thin1000.top0.01logp.subset.csv", row.names = FALSE, quote=FALSE)

########## PART 2: LFMM

# Select best run for K=2 from the snmf models
best <- which.min(cross.entropy(project, K=2))

# reading in sample info - maybe don't need this??
#sample_info <- read.csv("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current/narwhal/sample_info_narwhal_n57.csv", header=T)

# make lfmm file of the snps
ped2lfmm("narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.ped")

## Run lfmm! This is set up for doing 1 environmental variable at a time, then have to re-run this section for each variable

# manually made .env files for each environmental variable from the sample_info_narwhal_with_climdat_100KM.csv file

# Set up variables
lfmm_data <- "narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.lfmm"
lfmm_env <- "env files/narwhal_sst_100KM.env"

# Run lfmm2 with selected K
mod <- lfmm2(input=lfmm_data, env=lfmm_env, K=2)

# GEA significance test? showing K=2 estimated factors
#plot(mod@U, col = "grey", pch = 19, xlab = "Factor 1", ylab = "Factor 2")

# Computing P-values and plotting their minus log10 values 
pv <- lfmm2.test(object=mod, 
                 input=lfmm_data, 
                 env=lfmm_env,
                 linear=TRUE, genomic.control=TRUE)

# saving workspace so don't have to re-run later
save.image("narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000_sst_100KM.RData")

# re-load if needed if resuming from here
load("narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000_sst_100KM.RData")

# Let's work on adding snp info to the p values
snp_info <- read.table("narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.map")

# pull out only CHROM and POS
snp_info <- select(snp_info, V2)
colnames(snp_info) <- "snp"

pvals <- data.frame("snp"=snp_info,"zscore"=pv$zscores,"pvalue"=pv$pvalues)
# save this for whole snps set

# Use qvalue to use FDR rave 0.05 to mark significant snps
pvals$qval <- qvalue(pvals$pvalue, fdr.level=0.05)$signif
lfmm.snps <- pvals[pvals$qval == TRUE,"snp"]

# hardly any in this dataset, so let's save the whole thing instead
write.csv(pvals, "narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000_sst_100KM.pvalues.csv")

# manhattan plot
plot(-log10(pv$pvalues), col="grey", cex=.4, pch=19, main="sst_100KM_present")

# remove NAs?
#pvals <- pvals[! is.na(pvals$pvalue),]


########## PART 3:
## after running lfmm2 for each thing separately, next going to put results together into one dataframe.
# not many sig snps there either so let's go with the top 1% again for each
# then make list of unique snps (could be some overlap so want to avoid, OR maybe it doesn't matter if the list for filtering later has overlap?)

# moved stuff from desktop back into dropbox folders
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/gradient_forest_current/narwhal/LEA_narwhal_GF_rerun")

# read in top0.01 snps from genome scan:
scan <- read.csv("narwhal_K2.imputed.thin1000.top0.01logp.subset.csv")

# read in csv for lfmms (all snps, still need to filter top 1%) - code set up for doing this one at a time - maybe this part should've been set up in Part2
#sst <- read.csv("narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000_sst_100KM.pvalues.csv")
#icethick <- read.csv("narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000_icethick_100KM.pvalues.csv")
#salinity <- read.csv("narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000_salinity_100KM.pvalues.csv")
chloro <- read.csv("narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000_chloro_100KM.pvalues.csv")


# Not very efficient code but let's do the top0.01 extraction like we did for the snmf scan.
# one at a time, sst, ice thick, ...

# add snp info to the p values
snp_info <- read.table("narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.map")

# pull out only CHROM and POS
snp_info <- select(snp_info, V1, V4)

# get the log p's, add snp_info, order from large to small, pull top 1%
logp <- as.data.frame(-log10(chloro$pvalue))
logp$original_order <- 1:nrow(logp) #order number can also act as snp id later since this is for all 1,051,600 snps
snps_logp <- cbind(snp_info, logp)
colnames(snps_logp) <- c("CHROM", "POS", "logp", "order")
snps_logp2 <- snps_logp[,c(4,1,2,3)] #reordering the columns

# putting in order of logp (higher more significant)
snps_logp_order <- snps_logp2[order(-snps_logp2$logp),]

# set variable for number of snpts
count <- nrow(snps_logp_order)

# extract top x%, using 1% here
top_count <- floor(0.01*count)
snp_subset <- as.data.frame(snps_logp_order[1:top_count,])

topsnps <- snp_subset

# save it as csv too so dont have to redo later
write.csv(topsnps, "narwhal_K2.imputed.thin1000.lfmm.chloro_100KM.top0.01logp.subset.csv")

test <- snp_subset[order(snp_subset$order),]
plot((test$logp), col = "grey", cex = .4, pch = 19, main="chloro_100KM_present")

####
# ok, now reload the top 0.01 snp sets

# read in top0.01 snps from genome scan and envs
scan <- read.csv("narwhal_K2.imputed.thin1000.top0.01logp.subset.csv")
sst <- read.csv("narwhal_K2.imputed.thin1000.lfmm.sst_100KM.top0.01logp.subset.csv")
icethick <- read.csv("narwhal_K2.imputed.thin1000.lfmm.icethick_100KM.top0.01logp.subset.csv")
salinity <- read.csv("narwhal_K2.imputed.thin1000.lfmm.salinity_100KM.top0.01logp.subset.csv")
chloro <- read.csv("narwhal_K2.imputed.thin1000.lfmm.chloro_100KM.top0.01logp.subset.csv")
# current velocity doesn't really change in 2050 or 2100, maybe no point to add this in? exclude for now
#curvel <- read.csv("narwhal_K2.imputed.thin1000.lfmm.sst_100KM.top0.01logp.subset.csv")


#so 5 files with 10,523 snps each -> total of 52,615 SNPs

# let's see how many overlap? and also make a nice snp list to extract from vcf later

# just look at snp infos
library(tidyverse)

scan <- select(scan, -logp)
sst <- select(sst, c(-X,-logp))
icethick <- select(icethick, c(-X,-logp))
salinity <- select(salinity, c(-X,-logp))
chloro <- select(chloro, c(-X,-logp))

merged <- rbind(scan, sst, icethick, salinity, chloro)
merged_nodup <- merged[!duplicated(merged), ]

# ok so we have a list of unique snps for filtering vcf, 47,242 snps.
# we removed 5,373 snps that were in at least 2 things
# not sure yet how we want to look at overlaps
# anyway, continuing to map snp CHROM POS list for vcftools filtering

merged_nodup_order <- merged_nodup[order(merged_nodup$order),]

snp_list <- select(merged_nodup_order, -order)

write_delim(snp_list, file="narwhal_top0.01_scan_4env_merged_snplist.txt", delim = "\t")

