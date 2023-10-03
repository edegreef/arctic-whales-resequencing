# Continuing rehh stuff with narwhal and bowhead whale, using the output from data2haplohh step (02-haplohh_rehh_prep.R)

# Rehh vignette: https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html
# Help from Matt Thorstensen to set this up

# Load modules
library(rehh)
library(tidyverse)

# Working directory
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/bowhead/rehh")

# Load wgscan.rds file
wgscan <- readRDS("bowhead_wgscan.rds")

# Calculate genome-wide iHS values
wgscan.ihs <- ihh2ihs(wgscan)

# Default plots
freqbinplot(wgscan.ihs)
distribplot(wgscan.ihs$ihs$IHS, xlab = "iHS")
distribplot(wgscan.ihs$ihs$IHS, xlab = "iHS", qqplot = TRUE)
manhattanplot(wgscan.ihs, main = "iHS")
#manhattanplot(wgscan.ihs, pval = TRUE, threshold = 4, main = "p-value of iHS")

# Extracting ihs scores
ihs <- wgscan.ihs$ihs

# We want to use absolute values for ihs because we don't know which derived/ancestral status
ihs$IHS_ABS <- abs(ihs$IHS)

# Add snp ID
ihs$SNP <- paste(ihs$CHR, ihs$POSITION, sep="-")

#### For bowhead only: the scaffold automatically ordered like scaffold_0, scaffold_1, scaffold_11... want to just do scaffold0,1,2,3,4...
ihs[c("scaffold", "scafnum")] <- str_split_fixed(ihs$CHR, '_',2)
ihs$scafnumnum <- as.numeric(ihs$scafnum)
ihs <- ihs[order(ihs[,"scafnumnum"], ihs[,"POSITION"]),] # order by "scafnum" then by "POSITION"
####

# Back to pipeline for both species datasets
# Add plot_order row for plotting in ggplot
ihs$plot_order <- 1:nrow(ihs)


## Calculate candidate regions from the wgscan.ihs
ihs_candidates <- calc_candidate_regions(wgscan.ihs,
                                         threshold=4,
                                         pval=TRUE,
                                         window_size=1E3,#1E5 narwhal. #1E3 for bowhead
                                         overlap=1E2, #1E4 narwhal. 1E2 for bowhead
                                         min_n_extr_mrk=3)

ihs_candidates

# Filter for significant hits (q < 0.01) among the NOR_SS XP-EHH SNPs
# q < 0.05, 1.30103
sig_ihs <- filter(ihs, LOGPVALUE >= 4)

# Check which SNPs appear within a candidate region
sig_ihs$cand_region <- ifelse(sapply(sig_ihs$POSITION, function(p) 
  any(ihs_candidates$START <= p & ihs_candidates$END >= p)),"cand", NA)

ihs_cands <- subset(sig_ihs, cand_region=="cand")

# Save the ihs_cands
write.table(ihs_cands, "bowhead_ihs_cands_snps_nextrmark3_order.txt", sep="\t",quote=FALSE, row.names=FALSE)

# Extract 3 columns for SNPLOC file
snploc <- ihs_cands[,c("SNP", "CHR", "POSITION")]
write.table(snploc, "bowhead_ihs_cands_snps_nextrmark3_snploc.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)


rehh_plot <- ggplot()+
  geom_point(data=ihs, aes(x=plot_order,y=IHS_ABS,colour=CHR), alpha=0.5, size=1.5)+
  # Bowhead Whale
  scale_color_manual(values=rep(c("#dcb68f","gray60"), ceiling(length(ihs$CHR)/2))[1:length(ihs$CHR)])+
  geom_point(data=ihs_cands, aes(x=plot_order, y=IHS_ABS), alpha=0.8, size=2,colour="orange") +
  
  # Narwhal
  #scale_color_manual(values=rep(c("#a2c4d4","gray60"), ceiling(length(ihs$CHR)/2))[1:length(ihs$CHR)])+
  #geom_point(data=ihs_cands, aes(x=plot_order, y=IHS_ABS), alpha=0.8, size=2,colour="#2171B5") +
  
  # Other plot things
  theme_bw()+
  theme(legend.position="none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), 
        #panel.grid.minor.y=element_blank(),
        #panel.grid.major.y=element_blank(),
        text=element_text(size=14, family="serif"))+
       # theme(axis.text.x = element_text(angle = -45, hjust=1)))+
       # axis.text.x=element_blank())+
 # scale_x_reverse()+ # for narwhal, b/c for some reason it flipped order earlier
 # scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(limits=c(2,7.5))+ #(2,6.8) for narwhal, (2,7.5) for bowhead
  ylab("|iHS|")+
  xlab("Position")
 # geom_hline(yintercept=4, color="gray50", linetype="dashed")
  
rehh_plot

#ggsave("narwhal_IHS_ABS_pvalue4_plot2_mark3_blues_reverse.png", width=10, height=3, dpi=1000)
ggsave("bowhead_IHS_ABS_pvalue4_plot2_mark3_check_order.png", width=10, height=3, dpi=1000)

# Maybe save the ihs's too so easier to look at again later
write.csv(ihs, "bowhead_ihs_data.csv", row.names = F)

# If need to re-load
ihs <- read.csv("bowhead_ihs_data.csv")
