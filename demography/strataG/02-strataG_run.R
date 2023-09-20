# Estimating Ne with strataG

# install strataG
#devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)

# load libraries
library(vcfR)
library(adegenet)
library(strataG)

setwd("C:/Users/eveli/Dropbox/Whales with Garroway/strataG/narwhal")

# load vcf 
snpsR <- read.vcfR("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.n57.n20.nomiss.maf05.thinned25K.set1.vcf", verbose = T)

# convert vcf to genpop
snps_genind <- vcfR2genind(snpsR)
#class(snps_genind)

# convert genind to gtypes
snps_gtypes <- genind2gtypes(snps_genind)
#class(snps_gtypes)

# Estimating Ne using ldNe (from https://github.com/jdalapicolla/Ne_StrataG.R/blob/master/Ne_Estimation.R)
Ne <- ldNe(snps_gtypes, 
           maf.threshold=0, 
           by.strata=TRUE, 
           ci=0.95, 
           drop.missing=TRUE,
           num.cores=4)
Ne

# 25K snps over 16 samples took 4 minutes
# 50K snps over 16 samples start 11:23- after 20 minutes crashed/slowed down computer too much


# Save the results:
write.csv(Ne, "NeResults.narhwal.n57.25K.set1.csv")
