# Make a SFS file from VCF, using the frq.count file from vcftools

# Note for vcftools step is to run "vcftools --counts" on linux to get allele count info. 
# Also note, using vcf file that is filtered up to autosome step (NOT filtered for MAF and NOT LD pruned)
# Example: (here also removing the 1 kin pair individual in same line) 
# vcftools --gzvcf BOW_SNPS.filter1.miss.biallel.min100kb.autosomes.vcf.gz --remove-indv 99_01 --counts --out BOW_SNPS.filter1.miss.biallel.min100kb.autosomes.n22


library(readr)
library(tidyverse)
library(ggplot2)
library(reshape2)

# Load the frq.count output and rename columns
counts_raw <- read_delim("BOW_SNPS.filter1.miss.biallel.min100kb.autosomes.n22.frq.count", delim="\t",col_names=c("CHROM", "POS", "N_ALLELES", "N_CHR", "ALLELE_COUNT1","ALLELE_COUNT2"), skip = 1)

# First try with snps that do not have any missing data (N_CHR=44 in this case b/c 22 BOW samples and 2 alleles per site)
counts <- subset(counts_raw, N_CHR == 44)

# Need to edit the ALLELE_COUNT columns to extract the counts

# Remove the first two letters for all values in ALLELE_COUNT1 and ALLELE_COUNT2
# (probably more elegant ways to do this, but here using split function with ":")
split <- str_split_fixed(counts$ALLELE_COUNT1, ":", 2)
counts$ALLELE_COUNT1_NUM <- split[,2]

split2 <- str_split_fixed(counts$ALLELE_COUNT2, ":", 2)
counts$ALLELE_COUNT2_NUM <- split2[,2]

# Update dataframe with just the info we need
counts <- counts[,c("CHROM", "POS", "ALLELE_COUNT1_NUM", "ALLELE_COUNT2_NUM")]

# Data will be folded because ancestral allele info is unkown, thus we extract lower value count for each snp. Make another column with the lower value between ALLELE_COUNT1_NUM and ALLELE_COUNT2_NUM for each row

# Need to make them numeric first
counts$ALLELE_COUNT1_NUM <- as.numeric(counts$ALLELE_COUNT1_NUM)
counts$ALLELE_COUNT2_NUM <- as.numeric(counts$ALLELE_COUNT2_NUM)
           
# Extract minimum/lower count
counts$min <- pmin(counts$ALLELE_COUNT1_NUM, counts$ALLELE_COUNT2_NUM)

# Count number of occurrences in the min column
sfs <- as.data.frame(table(counts$min))

# Save into as SFS file
write.table(sfs, file="BOW_folded_SFS.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Plot folded sfs as histogram/bar-plot
ggplot(sfs, aes(y=Freq, x=Var1)) + 
  geom_bar(position="dodge", width=0.7, stat="identity",)+
  theme_classic()+
  scale_fill_manual(values=c("darkblue", "gray"))+
  xlab("allele frequency")+
  ylab("snps")

ggsave("BOW_folded_SFS_bar.png", width=5, height=3, dpi=200)
