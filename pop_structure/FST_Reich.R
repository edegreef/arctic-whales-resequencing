# Pairwise FST using Reich's unbiased estimate to handle unequal and small sample sizes. (from Reich et al. 2009: https://www.nature.com/articles/nature08365)

# 1) Function code to create Reich.Fst from Gaetano et al. 2014: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0091237#s6, supplement: https://doi.org/10.1371/journal.pone.0091237.s002)

Reich.Fst <- function(pop1,pop2,call.rate = 0.95, top.number = 10){
  # remove the SNPs that are not in common between the 2 populations
  snp.to.keep <- intersect(row.names(pop1),row.names(pop2))
  if (length(snp.to.keep) == 0){print("Error: no SNP in common");return(NULL)}
  pop1.k <- pop1[snp.to.keep,]
  pop2.k <- pop2[snp.to.keep,]
  # change the reference allele if is not concordant between the 2 populations
  if (sum(pop1.k$A1 == pop2.k$A1) != length(snp.to.keep)){
    idx <- which(pop1.k$A1 != pop2.k$A1)
    idx.rev <- which(pop1.k$A1 != pop2.k$A1 & pop1.k$A1 == pop2.k$A2)
    idx.rm  <- which(pop1.k$A1 != pop2.k$A1 & pop1.k$A1 != pop2.k$A2)
    if(length(idx.rev) > 0){
      provv <- pop1.k$A1[idx.rev]
      pop1.k$A1[idx.rev] <- pop1.k$A2[idx.rev]
      pop1.k$A2[idx.rev] <- provv
      provv <- pop1.k$x0[idx.rev]
      pop1.k$x0[idx.rev] <- pop1.k$x2[idx.rev]
      pop1.k$x2[idx.rev] <- provv}
    if(length(idx.rm) > 0){      
      pop1.k <- pop1.k[-idx.rm,]
      pop2.k <- pop2.k[-idx.rm,]}}
  # remove SNPs with low call rate in one or both populations
  x0 <- pop1.k$x0
  x1 <- pop1.k$x1
  x2 <- pop1.k$x2
  s <- x0 + x1 + x2
  y0 <- pop2.k$x0
  y1 <- pop2.k$x1
  y2 <- pop2.k$x2
  t <- y0 + y1 + y2
  idx.rm.pop1 <- which(s < max(s)*call.rate)
  idx.rm.pop2 <- which(t < max(t)*call.rate)
  idx.rm.all <- union(idx.rm.pop1,idx.rm.pop2)
  x0 <- x0[-idx.rm.all]
  x1 <- x1[-idx.rm.all]
  x2 <- x2[-idx.rm.all]
  s  <- s[-idx.rm.all]
  y0 <- y0[-idx.rm.all]
  y1 <- y1[-idx.rm.all]
  y2 <- y2[-idx.rm.all]
  t  <- t[-idx.rm.all]
  # compute SNP_Fst and global Fst estimators in presence of inbreeding   
  e.x <- ((x1 + 2*x2)/(2*s) - (y1 + 2*y2)/(2*t))^2 + x1/(4*s*s) + y1/(4*t*t)
  e.h1 <- (x0*x2 + (x0 + x2)*x1/2 + x1*(x1-1)/4)/(s*(s-1))
  e.h2 <- (y0*y2 + (y0 + y2)*y1/2 + y1*(y1-1)/4)/(t*(t-1))
  N <- e.x - e.h1/s - e.h2/t
  D <- N + e.h1 + e.h2
  Fst.v <- N/D
  names(Fst.v) <- row.names(pop1.k[-idx.rm.all,])
  Fst.o <- Fst.v[order(Fst.v,decreasing=TRUE)]
  F.global <- sum(N)/sum(D)
  se1 <- sd(N)/sqrt(length(N))
  se2 <- sd(D)/sqrt(length(N))
  se.F <- sqrt(se1*se1 + se2*se2)
  F_L95 <- F.global - 1.96*se.F
  F_U95 <- F.global + 1.96*se.F
  Z <- F.global/se.F
  p <- 2*(1 - pnorm(Z))
  if(p < 2e-16){p <- "less than 2e-16"}
  output <- list()
  output[[1]] <- c(F.global,F_L95,F_U95,p)
  names(output[[1]]) <- c("Reich.Fst","L.95%.CI","U.95%.CI","p.val")
  output[[2]] <- data.frame(Fst.o[1:top.number])
  names(output[[2]]) <- c("Reich.Fst")
  return(output)}


# 2) Set up data to run Reich's fst
# I couldn't get it to work as a combined genind object, so I did each pair separately. Formatting data into pop matrices with ncessary input info.

# From the supplemental page in Gaetano et al. 2014:
# "input data frame pop1 is a N x 5 matrix
# where N is the number of SNPs
# row names correspond to the SNP name
# x0 represent the number of samples with 0 copies of the variant allele
# x1 represent the number of samples with 1 copy of the variant allele
# x2 represent the number of samples with 2 copies of the variant allele
# A1 common allele
# A2 variant allele"


library(adegenet)
library(vcfR)
library(tidyverse)

setwd("C:/Users/eveli/Dropbox/Whales with Garroway/narwhal/snps")

# Read in vcfs for each pop
# Dropping IG, PB, and PG b/c they're 1-2 samples each
vcf_AB <- read.vcfR("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.AB.vcf.gz") 
vcf_BI <- read.vcfR("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.BI.vcf.gz") 
vcf_CR <- read.vcfR("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.CR.vcf.gz") 
vcf_GF <- read.vcfR("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.GF.vcf.gz") 
#vcf_IG <- read.vcfR("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.IG.vcf.gz") 
vcf_PB <- read.vcfR("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.PB.vcf.gz") 
#vcf_PG <- read.vcfR("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.PG.vcf.gz") 
vcf_PI <- read.vcfR("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.PI.vcf.gz") 
vcf_RB <- read.vcfR("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.RB.vcf.gz") 
vcf_RE <- read.vcfR("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.RE.vcf.gz") 
vcf_SB <- read.vcfR("narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.SB.vcf.gz") 


# Convert vcf to genlight object
genlight_AB <- vcfR2genlight(vcf_AB)
genlight_BI <- vcfR2genlight(vcf_BI)
genlight_CR <- vcfR2genlight(vcf_CR)
genlight_GF <- vcfR2genlight(vcf_GF)
genlight_IG <- vcfR2genlight(vcf_IG)
genlight_PB <- vcfR2genlight(vcf_PB)
genlight_PG <- vcfR2genlight(vcf_PG)
genlight_PI <- vcfR2genlight(vcf_PI)
genlight_RB <- vcfR2genlight(vcf_RB)
genlight_RE <- vcfR2genlight(vcf_RE)
genlight_SB <- vcfR2genlight(vcf_SB)

# Reformat data to matrices
mat_AB <- as.matrix(genlight_AB)
mat_BI <- as.matrix(genlight_BI)
mat_CR <- as.matrix(genlight_CR)
mat_GF <- as.matrix(genlight_GF)
mat_IG <- as.matrix(genlight_IG)
mat_PB <- as.matrix(genlight_PB)
mat_PG <- as.matrix(genlight_PG)
mat_PI <- as.matrix(genlight_PI)
mat_RB <- as.matrix(genlight_RB)
mat_RE <- as.matrix(genlight_RE)
mat_SB <- as.matrix(genlight_SB)


# Transpose matrix so rows=snps, and also convert to dataframe
mat_AB <- as.data.frame(t(mat_AB))
mat_BI <- as.data.frame(t(mat_BI))
mat_CR <- as.data.frame(t(mat_CR))
mat_GF <- as.data.frame(t(mat_GF))
mat_IG <- as.data.frame(t(mat_IG))
mat_PB <- as.data.frame(t(mat_PB))
mat_PG <- as.data.frame(t(mat_PG))
mat_PI <- as.data.frame(t(mat_PI))
mat_RB <- as.data.frame(t(mat_RB))
mat_RE <- as.data.frame(t(mat_RE))
mat_SB <- as.data.frame(t(mat_SB))


# Make columns of allele counts. The n in 'rowSums(mat[,1:n]...)' is sample size. (i.e. arctic has 11 samples). If the rowSums is not capped to sample size then it may include x0 - x2 columns in counts which will mess up counts.
# do for each pop
mat_AB$x0 <- rowSums(mat_AB[,1:9] == 0, na.rm=TRUE)
mat_AB$x1 <- rowSums(mat_AB[,1:9] == 1, na.rm=TRUE)
mat_AB$x2 <- rowSums(mat_AB[,1:9] == 2, na.rm=TRUE)

mat_BI$x0 <- rowSums(mat_BI[,1:4] == 0, na.rm=TRUE)
mat_BI$x1 <- rowSums(mat_BI[,1:4] == 1, na.rm=TRUE)
mat_BI$x2 <- rowSums(mat_BI[,1:4] == 2, na.rm=TRUE)

mat_CR$x0 <- rowSums(mat_CR[,1:3] == 0, na.rm=TRUE)
mat_CR$x1 <- rowSums(mat_CR[,1:3] == 1, na.rm=TRUE)
mat_CR$x2 <- rowSums(mat_CR[,1:3] == 2, na.rm=TRUE)

mat_GF$x0 <- rowSums(mat_GF[,1:10] == 0, na.rm=TRUE)
mat_GF$x1 <- rowSums(mat_GF[,1:10] == 1, na.rm=TRUE)
mat_GF$x2 <- rowSums(mat_GF[,1:10] == 2, na.rm=TRUE)

mat_IG$x0 <- rowSums(mat_IG[,1:2] == 0, na.rm=TRUE)
mat_IG$x1 <- rowSums(mat_IG[,1:2] == 1, na.rm=TRUE)
mat_IG$x2 <- rowSums(mat_IG[,1:2] == 2, na.rm=TRUE)

mat_PB$x0 <- rowSums(mat_PB[,1:2] == 0, na.rm=TRUE)
mat_PB$x1 <- rowSums(mat_PB[,1:2] == 1, na.rm=TRUE)
mat_PB$x2 <- rowSums(mat_PB[,1:2] == 2, na.rm=TRUE)

mat_PG$x0 <- rowSums(mat_PG[,1:1] == 0, na.rm=TRUE)
mat_PG$x1 <- rowSums(mat_PG[,1:1] == 1, na.rm=TRUE)
mat_PG$x2 <- rowSums(mat_PG[,1:1] == 2, na.rm=TRUE)

mat_PI$x0 <- rowSums(mat_PI[,1:4] == 0, na.rm=TRUE)
mat_PI$x1 <- rowSums(mat_PI[,1:4] == 1, na.rm=TRUE)
mat_PI$x2 <- rowSums(mat_PI[,1:4] == 2, na.rm=TRUE)

mat_RB$x0 <- rowSums(mat_RB[,1:9] == 0, na.rm=TRUE)
mat_RB$x1 <- rowSums(mat_RB[,1:9] == 1, na.rm=TRUE)
mat_RB$x2 <- rowSums(mat_RB[,1:9] == 2, na.rm=TRUE)

mat_RE$x0 <- rowSums(mat_RE[,1:9] == 0, na.rm=TRUE)
mat_RE$x1 <- rowSums(mat_RE[,1:9] == 1, na.rm=TRUE)
mat_RE$x2 <- rowSums(mat_RE[,1:9] == 2, na.rm=TRUE)

mat_SB$x0 <- rowSums(mat_SB[,1:4] == 0, na.rm=TRUE)
mat_SB$x1 <- rowSums(mat_SB[,1:4] == 1, na.rm=TRUE)
mat_SB$x2 <- rowSums(mat_SB[,1:4] == 2, na.rm=TRUE)

# Subset matrix to just get snp name and allele counts
allele_count_AB <- mat_AB %>% select(c(x0,x1,x2))
allele_count_BI <- mat_BI %>% select(c(x0,x1,x2))
allele_count_CR <- mat_CR %>% select(c(x0,x1,x2))
allele_count_GF <- mat_GF %>% select(c(x0,x1,x2))
allele_count_IG <- mat_IG %>% select(c(x0,x1,x2))
allele_count_PB <- mat_PB %>% select(c(x0,x1,x2))
#allele_count_PG <- mat_PG %>% select(c(x0,x1,x2)) # dont think i can use PG with only n=1
allele_count_PI <- mat_PI %>% select(c(x0,x1,x2))
allele_count_RB <- mat_RB %>% select(c(x0,x1,x2))
allele_count_RE <- mat_RE %>% select(c(x0,x1,x2))
allele_count_SB <- mat_SB %>% select(c(x0,x1,x2))

# Get the allele infos (A/T/G/C), read vcfs as tables
table_AB <- read.table("NAR_SNPS.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.AB.vcf.gz")
table_BI <- read.table("NAR_SNPS.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.BI.vcf.gz")
table_CR <- read.table("NAR_SNPS.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.CR.vcf.gz")
table_GF <- read.table("NAR_SNPS.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.GF.vcf.gz")
table_IG <- read.table("NAR_SNPS.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.IG.vcf.gz")
table_PB <- read.table("NAR_SNPS.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.PB.vcf.gz")
table_PI <- read.table("NAR_SNPS.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.PI.vcf.gz")
table_RB <- read.table("NAR_SNPS.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.RB.vcf.gz")
table_RE <- read.table("NAR_SNPS.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.RE.vcf.gz")
table_SB <- read.table("NAR_SNPS.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n57.SB.vcf.gz")



# Extract the common allele (A1) and variant allele (A2) 
variants_AB <- table_AB %>% select(c(V4,V5))
colnames(variants_AB) <- c("A1", "A2")

variants_BI <- table_BI %>% select(c(V4,V5))
colnames(variants_BI) <- c("A1", "A2")

variants_CR <- table_CR %>% select(c(V4,V5))
colnames(variants_CR) <- c("A1", "A2")

variants_GF <- table_GF %>% select(c(V4,V5))
colnames(variants_GF) <- c("A1", "A2")

variants_IG <- table_IG %>% select(c(V4,V5))
colnames(variants_IG) <- c("A1", "A2")

variants_PB <- table_PB %>% select(c(V4,V5))
colnames(variants_PB) <- c("A1", "A2")

variants_PI <- table_PI %>% select(c(V4,V5))
colnames(variants_PI) <- c("A1", "A2")

variants_RB <- table_RB %>% select(c(V4,V5))
colnames(variants_RB) <- c("A1", "A2")

variants_RE <- table_RE %>% select(c(V4,V5))
colnames(variants_RE) <- c("A1", "A2")

variants_SB <- table_SB %>% select(c(V4,V5))
colnames(variants_SB) <- c("A1", "A2")


# Merge allele counts and types, then it should be in the right format for Reich.Fst 
pop_AB <- cbind(allele_count_AB, variants_AB)
pop_BI <- cbind(allele_count_BI, variants_BI)
pop_CR <- cbind(allele_count_CR, variants_CR)
pop_GF <- cbind(allele_count_GF, variants_GF)
pop_IG <- cbind(allele_count_IG, variants_IG)
pop_PB <- cbind(allele_count_PB, variants_PB)
pop_PI <- cbind(allele_count_PI, variants_PI)
pop_RB <- cbind(allele_count_RB, variants_RB)
pop_RE <- cbind(allele_count_RE, variants_RE)
pop_SB <- cbind(allele_count_SB, variants_SB)



# Run Reich's fst for each pair
# example:
#fst_AR_IC <- Reich.Fst(pop_arctic, pop_iceland, call.rate=0.75, top.number=10)

# AB, BI, CR, GF, PI, RB, RE, SB
fst_AB_BI <- Reich.Fst(pop_AB, pop_BI, call.rate=0.75, top.number=10)
fst_AB_CR <- Reich.Fst(pop_AB, pop_CR, call.rate=0.75, top.number=10)
fst_AB_GF <- Reich.Fst(pop_AB, pop_GF, call.rate=0.75, top.number=10)
fst_AB_IG <- Reich.Fst(pop_AB, pop_IG, call.rate=0.75, top.number=10)
fst_AB_PB <- Reich.Fst(pop_AB, pop_PB, call.rate=0.75, top.number=10)
fst_AB_PI <- Reich.Fst(pop_AB, pop_PI, call.rate=0.75, top.number=10)
fst_AB_RB <- Reich.Fst(pop_AB, pop_RB, call.rate=0.75, top.number=10)
fst_AB_RE <- Reich.Fst(pop_AB, pop_RE, call.rate=0.75, top.number=10)
fst_AB_SB <- Reich.Fst(pop_AB, pop_SB, call.rate=0.75, top.number=10)

fst_BI_CR <- Reich.Fst(pop_BI, pop_CR, call.rate=0.75, top.number=10)
fst_BI_GF <- Reich.Fst(pop_BI, pop_GF, call.rate=0.75, top.number=10)
fst_BI_IG <- Reich.Fst(pop_BI, pop_IG, call.rate=0.75, top.number=10)
fst_BI_PB <- Reich.Fst(pop_BI, pop_PB, call.rate=0.75, top.number=10)
fst_BI_PI <- Reich.Fst(pop_BI, pop_PI, call.rate=0.75, top.number=10)
fst_BI_RB <- Reich.Fst(pop_BI, pop_RB, call.rate=0.75, top.number=10)
fst_BI_RE <- Reich.Fst(pop_BI, pop_RE, call.rate=0.75, top.number=10)
fst_BI_SB <- Reich.Fst(pop_BI, pop_SB, call.rate=0.75, top.number=10)

fst_CR_GF <- Reich.Fst(pop_CR, pop_GF, call.rate=0.75, top.number=10)
fst_CR_IG <- Reich.Fst(pop_CR, pop_IG, call.rate=0.75, top.number=10)
fst_CR_PB <- Reich.Fst(pop_CR, pop_PB, call.rate=0.75, top.number=10)
fst_CR_PI <- Reich.Fst(pop_CR, pop_PI, call.rate=0.75, top.number=10)
fst_CR_RB <- Reich.Fst(pop_CR, pop_RB, call.rate=0.75, top.number=10)
fst_CR_RE <- Reich.Fst(pop_CR, pop_RE, call.rate=0.75, top.number=10)
fst_CR_SB <- Reich.Fst(pop_CR, pop_SB, call.rate=0.75, top.number=10)

fst_GF_IG <- Reich.Fst(pop_GF, pop_IG, call.rate=0.75, top.number=10)
fst_GF_PB <- Reich.Fst(pop_GF, pop_PB, call.rate=0.75, top.number=10)
fst_GF_PI <- Reich.Fst(pop_GF, pop_PI, call.rate=0.75, top.number=10)
fst_GF_RB <- Reich.Fst(pop_GF, pop_RB, call.rate=0.75, top.number=10)
fst_GF_RE <- Reich.Fst(pop_GF, pop_RE, call.rate=0.75, top.number=10)
fst_GF_SB <- Reich.Fst(pop_GF, pop_SB, call.rate=0.75, top.number=10)

fst_IG_PB <- Reich.Fst(pop_IG, pop_PB, call.rate=0.75, top.number=10)
fst_IG_PI <- Reich.Fst(pop_IG, pop_PI, call.rate=0.75, top.number=10)
fst_IG_RB <- Reich.Fst(pop_IG, pop_RB, call.rate=0.75, top.number=10)
fst_IG_RE <- Reich.Fst(pop_IG, pop_RE, call.rate=0.75, top.number=10)
fst_IG_SB <- Reich.Fst(pop_IG, pop_SB, call.rate=0.75, top.number=10)

fst_PB_PI <- Reich.Fst(pop_PB, pop_PI, call.rate=0.75, top.number=10)
fst_PB_RB <- Reich.Fst(pop_PB, pop_RB, call.rate=0.75, top.number=10)
fst_PB_RE <- Reich.Fst(pop_PB, pop_RE, call.rate=0.75, top.number=10)
fst_PB_SB <- Reich.Fst(pop_PB, pop_SB, call.rate=0.75, top.number=10)

fst_PI_RB <- Reich.Fst(pop_PI, pop_RB, call.rate=0.75, top.number=10)
fst_PI_RE <- Reich.Fst(pop_PI, pop_RE, call.rate=0.75, top.number=10)
fst_PI_SB <- Reich.Fst(pop_PI, pop_SB, call.rate=0.75, top.number=10)

fst_RB_RE <- Reich.Fst(pop_RB, pop_RE, call.rate=0.75, top.number=10)
fst_RB_SB <- Reich.Fst(pop_RB, pop_SB, call.rate=0.75, top.number=10)

fst_RE_SB <- Reich.Fst(pop_RE, pop_SB, call.rate=0.75, top.number=10)


# Save results in a thing
a <- data.frame(matrix(unlist(fst_AB_BI[1]), nrow=length(fst_AB_BI[1]), byrow=TRUE))
b <- data.frame(matrix(unlist(fst_AB_CR[1]), nrow=length(fst_AB_CR[1]), byrow=TRUE))
c <- data.frame(matrix(unlist(fst_AB_GF[1]), nrow=length(fst_AB_GF[1]), byrow=TRUE))
d <- data.frame(matrix(unlist(fst_AB_IG[1]), nrow=length(fst_AB_IG[1]), byrow=TRUE))
e <- data.frame(matrix(unlist(fst_AB_PB[1]), nrow=length(fst_AB_PB[1]), byrow=TRUE))
f <- data.frame(matrix(unlist(fst_AB_PI[1]), nrow=length(fst_AB_PI[1]), byrow=TRUE))
g <- data.frame(matrix(unlist(fst_AB_RB[1]), nrow=length(fst_AB_RB[1]), byrow=TRUE))
h <- data.frame(matrix(unlist(fst_AB_RE[1]), nrow=length(fst_AB_RE[1]), byrow=TRUE))
i <- data.frame(matrix(unlist(fst_AB_SB[1]), nrow=length(fst_AB_SB[1]), byrow=TRUE))

j <- data.frame(matrix(unlist(fst_BI_CR[1]), nrow=length(fst_BI_CR[1]), byrow=TRUE))
k <- data.frame(matrix(unlist(fst_BI_GF[1]), nrow=length(fst_BI_GF[1]), byrow=TRUE))
l <- data.frame(matrix(unlist(fst_BI_IG[1]), nrow=length(fst_BI_IG[1]), byrow=TRUE))
m <- data.frame(matrix(unlist(fst_BI_PB[1]), nrow=length(fst_BI_PB[1]), byrow=TRUE))
n <- data.frame(matrix(unlist(fst_BI_PI[1]), nrow=length(fst_BI_PI[1]), byrow=TRUE))
o <- data.frame(matrix(unlist(fst_BI_RB[1]), nrow=length(fst_BI_RB[1]), byrow=TRUE))
p <- data.frame(matrix(unlist(fst_BI_RE[1]), nrow=length(fst_BI_RE[1]), byrow=TRUE))
q <- data.frame(matrix(unlist(fst_BI_SB[1]), nrow=length(fst_BI_SB[1]), byrow=TRUE))

r <- data.frame(matrix(unlist(fst_CR_GF[1]), nrow=length(fst_CR_GF[1]), byrow=TRUE))
s <- data.frame(matrix(unlist(fst_CR_IG[1]), nrow=length(fst_CR_IG[1]), byrow=TRUE))
t <- data.frame(matrix(unlist(fst_CR_PB[1]), nrow=length(fst_CR_PB[1]), byrow=TRUE))
u <- data.frame(matrix(unlist(fst_CR_PI[1]), nrow=length(fst_CR_PI[1]), byrow=TRUE))
v <- data.frame(matrix(unlist(fst_CR_RB[1]), nrow=length(fst_CR_RB[1]), byrow=TRUE))
w <- data.frame(matrix(unlist(fst_CR_RE[1]), nrow=length(fst_CR_RE[1]), byrow=TRUE))
x <- data.frame(matrix(unlist(fst_CR_SB[1]), nrow=length(fst_CR_SB[1]), byrow=TRUE))

y <- data.frame(matrix(unlist(fst_GF_IG[1]), nrow=length(fst_GF_IG[1]), byrow=TRUE))
z <- data.frame(matrix(unlist(fst_GF_PB[1]), nrow=length(fst_GF_PB[1]), byrow=TRUE))
aa <- data.frame(matrix(unlist(fst_GF_PI[1]), nrow=length(fst_GF_PI[1]), byrow=TRUE))
ab <- data.frame(matrix(unlist(fst_GF_RB[1]), nrow=length(fst_GF_RB[1]), byrow=TRUE))
ac <- data.frame(matrix(unlist(fst_GF_RE[1]), nrow=length(fst_GF_RE[1]), byrow=TRUE))
ad <- data.frame(matrix(unlist(fst_GF_SB[1]), nrow=length(fst_GF_SB[1]), byrow=TRUE))

ae <- data.frame(matrix(unlist(fst_IG_PB[1]), nrow=length(fst_IG_PB[1]), byrow=TRUE))
af <- data.frame(matrix(unlist(fst_IG_PI[1]), nrow=length(fst_IG_PI[1]), byrow=TRUE))
ag <- data.frame(matrix(unlist(fst_IG_RB[1]), nrow=length(fst_IG_RB[1]), byrow=TRUE))
ah <- data.frame(matrix(unlist(fst_IG_RE[1]), nrow=length(fst_IG_RE[1]), byrow=TRUE))
ai <- data.frame(matrix(unlist(fst_IG_SB[1]), nrow=length(fst_IG_SB[1]), byrow=TRUE))

aj <- data.frame(matrix(unlist(fst_PB_PI[1]), nrow=length(fst_PB_PI[1]), byrow=TRUE))
ak <- data.frame(matrix(unlist(fst_PB_RB[1]), nrow=length(fst_PB_RB[1]), byrow=TRUE))
al <- data.frame(matrix(unlist(fst_PB_RE[1]), nrow=length(fst_PB_RE[1]), byrow=TRUE))
am <- data.frame(matrix(unlist(fst_PB_SB[1]), nrow=length(fst_PB_SB[1]), byrow=TRUE))

an <- data.frame(matrix(unlist(fst_PI_RB[1]), nrow=length(fst_PI_RB[1]), byrow=TRUE))
ao <- data.frame(matrix(unlist(fst_PI_RE[1]), nrow=length(fst_PI_RE[1]), byrow=TRUE))
ap <- data.frame(matrix(unlist(fst_PI_SB[1]), nrow=length(fst_PI_SB[1]), byrow=TRUE))

aq <- data.frame(matrix(unlist(fst_RB_RE[1]), nrow=length(fst_RB_RE[1]), byrow=TRUE))
ar <- data.frame(matrix(unlist(fst_RB_SB[1]), nrow=length(fst_RB_SB[1]), byrow=TRUE))

as <- data.frame(matrix(unlist(fst_RE_SB[1]), nrow=length(fst_RE_SB[1]), byrow=TRUE))

merged <- rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,ab,ac,ad,ae,af,ag,ah,ai,aj,ak,al,am,an,ao,ap,aq,ar,as)

# AB, BI, CR, GF, PI, RB, RE, SB
# making column with pair labels to add to merged data (make sure in same order)
pair <- data.frame(c(
  "AB-BI", "AB-CR", "AB-GF", "AB-IG", "AB-PB", "AB-PI", "AB-RB", "AB-RE", "AB-SB", 
  "BI-CR", "BI-GF", "BI-IG", "BI-PB",  "BI-PI", "BI-RB", "BI-RE", "BI-SB", 
  "CR-GF", "CR-IG", "CR-PB",  "CR-PI", "CR-RB", "CR-RE", "CR-SB", 
  "GF-IG", "GF-PB", "GF-PI", "GF-RB", "GF-RE", "GF-SB", 
  "IG-PB", "IG-PI", "IG-RB", "IG-RE", "IG-SB",
  "PB-PI", "PB-RB", "PB-RE", "PB-SB",
  "PI-RB", "PI-RE", "PI-SB", 
  "RB-RE", "RB-SB", 
  "RE-SB"))

fst_together <- cbind(pair, merged)
colnames(fst_together) <- c("pair", "Reich_fst", "ll_95CI", "ul_95CI", "p_val")

write.csv(fst_together, "Reichs_fst_output_NAR_10sites.csv")
