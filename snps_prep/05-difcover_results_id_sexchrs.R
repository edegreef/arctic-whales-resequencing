# Looking at DifCover output (.DNAcopyout) - https://github.com/timnat/DifCover
# Using male and female samples to determine X and Y linked scaffolds. Doing this for 4 male v female runs (M1F1, M1F2, M2F1, M2F2) and 2 controls (M1M2, F1F2). Script for the 2 controls towards the bottom.
# Guidance and code help from Phil Grayson (https://sexfindr.readthedocs.io/en/latest/) with some modifications

library(tidyverse)

setwd("C:/Users/Evelien de Greef/Dropbox/Whales with Garroway/BOW/difcover/control/M1F1")

# Enter in file infos
DNAcopyout ="M1_F1.ratio_per_w_CC0_a5_A45_b5_B48_v1000_l500.log2adj_1.067.DNAcopyout"
scaffold_info="C:/Users/Evelien de Greef/Dropbox/Whales with Garroway/BOW/difcover/BOW_reference.fasta.fai"

# Some more labels for later (plot title/label)
sample1="M1"
sample2="F1"
type="M1F1"

# Load in data, rename columns, and add bp spanned. Enrichment score: log2(sample1 coverage/sample2 coverage)
difcover <- read_tsv(file = DNAcopyout ,col_names = F) %>% 
  rename(scaf = X1, start = X2, stop = X3, windows = X4, enrichment_score = X5) %>% 
  mutate("bases spanned" = stop-start)

# Use fasta.fai for scaffold name and length
scaffold_lengths <- read_tsv(scaffold_info,col_names = c("scaf","length"))
scaffold_lengths <- scaffold_lengths[c("scaf", "length")]

# Join with difcover output
proportion <- full_join(difcover,scaffold_lengths) %>% 
  mutate(proportion = `bases spanned`/length)

# Initial plot of proportion of scaffold versus log2(male/female) coverage
proportion %>% 
  ggplot(aes(x=enrichment_score, y=proportion)) + 
  geom_point()+
  xlab("enrichment score")
  #geom_vline(xintercept=0, col="red")

# If want to save prelim plot:
#ggsave(paste("coverageplot_quick_", sample1, "_", sample2, "_", type, ".png",sep=""), width = 10, height = 6, dpi = 300)

# Interested in difCover regions that are likely to be Y (enrichment greater than 2) or X (enrichment less than -0.7369656)
filtered_proportion <- proportion %>% 
  filter(enrichment_score >= 2 | enrichment_score <= -0.7369656) %>% 
  group_by(scaf) %>%
  mutate("total scaf proportion with sig diff cov" = sum(`bases spanned`)/length)

# plot proportions to see how many are substantial proportions of their scaffold; based on plot, may want to only consider scaffolds that fit a cutoff like 0.25 (set as red line here)
filtered_proportion %>% 
  ggplot(aes(x=scaf,y=`total scaf proportion with sig diff cov`)) + 
  geom_point(alpha=0.2) + 
  geom_hline(yintercept=0.25, colour="red") + 
  ylab("total proportion") + 
  xlab("scaffold")+
  theme_classic()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

#ggsave(paste("scaffoldproportion_sig_", type, ".png",sep=""), width = 7, height = 5, dpi = 300)
# +geom_text(aes(label=scaf), alpha=0.5) 

# Keeping only scaffolds that have proportion of at least 0.25
filtered_proportion %>% 
  filter(`total scaf proportion with sig diff cov` >= 0.25) %>% 
  ggplot(aes(x=scaf,y=`total scaf proportion with sig diff cov`)) + 
  geom_point(alpha=0.2)+
  theme_classic()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

# number of scaffolds remaining using 0.25 cutoff
sex_linked <- filtered_proportion %>% 
  filter(`total scaf proportion with sig diff cov` >= 0.25)
nrow(sex_linked)

# Need to be careful about scaffolds that may have windows in X and other windows in Y.

# If many scaffolds, should check with # unique scaffolds instead of row numbers b/c of split windows
n_distinct(sex_linked$scaf)

# Number of scaffolds in X and Y
X_linked <- sex_linked %>% 
  filter(enrichment_score <= -0.7369656)
n_distinct(X_linked$scaf)

Y_linked <- sex_linked %>% 
  filter(enrichment_score >= 2) %>%
  filter(proportion >=0.25) # added this part because otherwise it may have a false hit if there was a scaffold on the X that had a short window that matched Y
n_distinct(Y_linked$scaf)

# Number of bases in X and Y
X_length <- sum(X_linked$`bases spanned`)
X_length

Y_length <- sum(Y_linked$`bases spanned`)
Y_length

# Proportion of X and Y in whole genome
genome_length <- sum(scaffold_lengths$length)
X_length/genome_length #something around 4-6% should be good for this species
Y_length/genome_length #Y is expected to be very small

# Extract scaffolds for X and Y into lists and give appropriate labels
X_linked_list <- X_linked %>% distinct(scaf)
X_linked_list$chr_type <- "X"

Y_linked_list <- Y_linked %>% distinct(scaf)
Y_linked_list$chr_type <- "Y"

linked_with_names <- bind_rows(X_linked_list,Y_linked_list)

# Combine that with the original "proportion" tibble and if there isn't a chr_type then label as autosome
proportion_with_names <- left_join(proportion,linked_with_names) %>% 
  replace_na(list(chr_type = "Autosome"))

# Plot proportion of scaffold versus log2(male/female) coverage with color codes
# Each dot is a window or set of windows in a scaffold
coverage_plot <- ggplot()+
  geom_point(data=proportion_with_names, aes(x=enrichment_score, y=proportion, color=chr_type))+
  xlab("enrichment score") + 
  ylab("proportion of scaffold")+
  theme_classic() + 
  scale_colour_manual(values=c("gray80","forestgreen","royalblue"))+
  ggtitle(paste("scaffold coverage b/t",sample1, "&", sample2, "-",type, sep=" "))#+ geom_vline(xintercept =0)

coverage_plot
ggsave(paste("coverageplot_", sample1, "_", sample2, "_", type, ".png",sep=""), width = 10, height = 6, dpi = 300)

# Prepare data for making bed files from the likely_X/Y_linked data (to use later for bedtools annotate)
X_linked$sexchr <- "X"
forbed_X <- X_linked %>% select(scaf, start, stop, sexchr)
write.table(forbed_X, file=paste("X_linked_scafwindows.", sample1, "&", sample2, "_", type, ".bed",sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

Y_linked$sexchr <- "Y"
forbed_Y <- Y_linked %>% select(scaf, start, stop, sexchr)
write.table(forbed_Y, paste("Y_linked_scafwindows.",sample1, "&", sample2, "_", type, ".bed",sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)



######## Looking at same-sex control data (M1M2, F1F2):

setwd("C:/Users/Evelien de Greef/Dropbox/Whales with Garroway/BOW/difcover/control/F1F2")

# Enter file infos
DNAcopyout ="BM_NSA_2009_BM_NSA_2012.ratio_per_w_CC0_a5_A48_b5_B42_v1000_l500.log2adj_0.875.DNAcopyout"
scaffold_info="C:/Users/Evelien de Greef/Dropbox/Whales with Garroway/BOW/difcover/BOW_reference.fasta.fai"

# Some more labels for later (plot title/label)
sample1="BM_NSA_2009_02"
sample2="BM_NSA_2012_02"
type="F1_F2"

# Load in the raw data (including non-significant regions), renaming columns, and adding bp spanned
difcover <- read_tsv(file = DNAcopyout ,col_names = F) %>% 
  rename(scaf=X1,start=X2,stop=X3,windows=X4,"enrichment_score"=X5) %>% 
  mutate("bases spanned" = stop-start)

# Parse samtools faidx output for scaffold name and length & join with difcover
scaffold_lengths <- read_tsv(scaffold_info,col_names = c("scaf","length"))
proportion <- full_join(difcover,scaffold_lengths) %>% 
  mutate(proportion = `bases spanned`/length)

# Plot proportion of scaffold versus log2(male/female) coverage
proportion %>% 
  ggplot(aes(x=enrichment_score,y=proportion)) + 
  geom_point()+xlab("enrichment_score")#+geom_vline(xintercept=0, col="white")
#ggsave(paste("coverage.", type, ".png",sep=""), width = 11, height = 8, dpi = 300)

# To keep not-enriched sites:
filtered_proportion_not_enrich <- proportion %>% 
  filter(enrichment_score < 2 & enrichment_score > -2) %>% 
  group_by(scaf) %>% 
  mutate("total scaf proportion with sig diff cov" = sum(`bases spanned`)/length)

# Plot not-enriched scafs
filtered_proportion_not_enrich %>% 
  ggplot(aes(x=scaf,y=`total scaf proportion with sig diff cov`)) + 
  geom_point(alpha=0.2) + 
  geom_hline(yintercept=0.25, colour="red") + 
  ggtitle(paste("total scaf proportion without sig different coverage -", type, sep=" ")) + 
  ylab("total proportion") + 
  xlab("scaffold")
#ggsave(paste("scaffoldproportion_notenrich_", type, ".png",sep=""), width = 7, height = 5, dpi = 300)

# Filter scaffolds with 0.25 cutoff
likely_not_enriched <- filtered_proportion_not_enrich %>% 
  filter(`total scaf proportion with sig diff cov` >= 0.25)

nrow(likely_not_enriched)
notenrich_length <- sum(likely_not_enriched$`bases spanned`)
notenrich_length
genome_length <- sum(scaffold_lengths$length)
notenrich_length/genome_length #should be high proportion, close to 100%

likely_not_enriched_list <- likely_not_enriched %>% 
  mutate(chr_type="Control") %>% 
  select(scaf,chr_type)

proportion_with_names <- left_join(proportion,likely_not_enriched_list) %>% 
  replace_na(list(chr_type = "Outlier"))

coverage_plot2 <- ggplot()+
  geom_point(data=proportion_with_names, aes(x=enrichment_score, y=proportion, color=chr_type))+
  xlab("enrichment_score") + 
  theme_classic() + 
  scale_colour_manual(values=c("gray80","tomato"))+
  ggtitle(paste("coverage plot of scaffolds between male and male sample -", type, sep=" ")) #+ geom_vline(xintercept =0)

coverage_plot2
ggsave(paste("coverageplot_", sample1, "_", sample2, "_", type, ".png",sep=""), width = 10, height = 6, dpi = 300)

# Prepare data for making bed file (to use later for bedtools annotate)
likely_not_enriched$chr <- "Control"
likely_not_enriched$rownum  <- 1:nrow(likely_not_enriched)
forbed_control <- likely_not_enriched %>% select(scaf, start, stop, chr, rownum)

# Before saving bed file for the control runs, check that no window "stop"s are at 0, otherwise bedtools will give error. Check the "likely_not_enriched" dataframe and sort the "stop" column in order by low-> high and then if there's a scaffold window with start and stop both at 0, then remove the row.

# F1F2: scaffold 132977 has a window 0 to 0 but also has 3 more windows in next line. Row 3515 in "forbed_control" should be removed (0 to 0). Not sure why this row exists, but I suspect this extra empty row happened in the difcover run.
# M1M2: scaffold 204716 same thing. Row 21806

#forbed_control <- forbed_control[-3515,-5] 
#forbed_control <- forbed_control[-21806,-5]
write.table(forbed_control, file=paste("Control_scafwindows",sample1, "&", sample2, "_", type, ".bed",sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Do same thing for M1M2 run

