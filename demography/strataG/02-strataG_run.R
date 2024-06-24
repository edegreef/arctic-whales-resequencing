# Estimating NE

#devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)

library(vcfR)
library(adegenet)
library(strataG)

setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/contemporary_Ne/strataG/narwhal")

# load vcf
snpsR <- read.vcfR("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.n57.10MB.nomiss.thinned25K.set1.vcf", verbose = T)

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

# Save the results:
write.csv(Ne, "NeResults.narwhal.n57.10MB.25K_nomaf.set1.csv")

