# Running PCA using help from the PCAdapt vignette: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html

# Load modules
library(pcadapt)
library(ggplot2)
library(patchwork)
library(dplyr)

# Working directory, go to directory for necessary data files
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/bowhead")

# Load data
snp_data <- read.pcadapt("rightwhale_map/bowhead_RWmap_snps_rmfix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n20.bed", type = "bed")

# Load sample info
sample_info <- read.csv("sample_info_bowhead.csv", header=T)

# Filter out individuals in the metadata that were removed in snp files
sample_info <- subset(sample_info, temp_remove!="X") 

# Verify snp file IDs matches order of metadata
# Load fam file
snp_IDs <- read.table("rightwhale_map/bowhead_RWmap_snps_rmfix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n20.fam")
snp_IDs$sample_ID <- sample_info$sample_ID
snp_IDs$all_equal <- snp_IDs$V1==snp_IDs$sample_ID #column should all say "TRUE"

# Run pcadapt. K value will be how many eigenvectors to be produced
x <- pcadapt(input = snp_data, K = 18)

# Screeplot
plot(x, option = "screeplot")

# Plot quick PCA
plot(x, option = "scores", pop = sample_info$location)

# Plot other eigenvectors
plot(x, option = "scores", i = 5, j = 6, pop = sample_info$location)

# Look at pca scores
scores <- as.data.frame(x$scores)

# Look at loadings
loadings <- as.data.frame(x$loadings)

# Z scores
z_scores <- as.data.frame(x$zscores)

# Look at proportion variance
proportion <- as.data.frame(x$singular.values)
proportion$squared <- proportion$`x$singular.values`* proportion$`x$singular.values`
prop_var <- as.data.frame(proportion$squared)
PC1_proportion <- (round(prop_var[1,], digits=4))*100
PC2_proportion <- (round(prop_var[2,], digits=4))*100
PC3_proportion <- (round(prop_var[3,], digits=4))*100
PC4_proportion <- (round(prop_var[4,], digits=4))*100
PC5_proportion <- (round(prop_var[5,], digits=4))*100
PC6_proportion <- (round(prop_var[6,], digits=4))*100

# Plot scree on own pref
prop_var$num <- 1:nrow(prop_var)
scree <- ggplot(data=prop_var, aes(x=num, y=prop_var$`proportion$squared`))+
  geom_point()+
  geom_line()+
  theme_bw()+
  ylab("Proportion of explained variance")+
  xlab("PC")
scree

#ggsave("scree_plot_LDprunedr08_n57.png", width=6, height=4.5, dpi=300)


#(nrow(x$scores) - 1) * length(x$pass)

# Save the scores as separate data file to adjust the pca plot for colors, etc
evec <- cbind(sample_info$sample_ID, scores)
colnames(evec)[1] <- "sample"
#write.csv(evec,"NBW2_SNPS_2M_pops_37_min50kb_LDprunedr08_standardID.evec.csv", row.names=FALSE)

library(ggrepel) # if want to label points to see which samples outliers

# Plot in ggplot
# Can switch V1 and V2 to different V's to look at other PCs- just remember to adjust the xlab and ylab
pca <- ggplot(data=evec, aes(x=V1,y=V2))+
  #previously used point size 3, but increasing 
  geom_point(aes(color=sample_info$location_ID, shape=sample_info$group, size=sample_info$group), alpha=0.8, size=3)+  
  theme_classic()+
# shapes for narwhal PCA
 # scale_shape_manual(values=c(17, 15, 18))+
 # scale_size_manual(values=c(3,3,4))+
  xlab(paste("PC1 (", PC1_proportion, "%)", sep=""))+
  ylab(paste("PC2 (", PC2_proportion, "%)", sep=""))+
  #geom_text_repel(aes(label=sample_info$sample_ID))+
  #geom_text_repel(label=sample_info$location_ID)+
  labs(color= "Region")+
  
  # # Narwhal colour scheme
 # scale_colour_manual(values=c("#54278F", "#807DBA", "#BCBDDC", "#0C2C84", "#225EA8", "#1D91C0", "#41B6C4", "#7FCDBB", "#FED976", "#FD8D3C", "#E31A1C"),             
 #                     breaks=c("GF", "RE", "SB", "IG", "PB", "AB", "PI", "CR", "BI", "PG", "RB"))
 # theme(text=element_text(family="serif"))

  # Bowhead colour scheme
 scale_color_manual(values=c("#807dbaff", "#e052c3ff","#ff9ae0ff", "#e31a1cff", "#965e35ff", "#7c0001ff", "#fd8d3cff","#fed976ff","#00b234ff" ,"#5ffeb6ff","#006531ff","#7fcdbbff","#41b6c4ff","#1d91c0ff", "#481edcff", "#0052ffff"),
         breaks=c("SB", "PB", "HB", "RB", "CH", "RI","CD", "WB", "IQ", "PG","NF","CR", "PI", "AB", "GH", "DIS"))

pca


#ggsave("pop_structure_results/PC1-2_bowhead_snps_LDprunedr08_n22_LABEL.png", width=10, height=8, dpi=1200)
#ggsave("pop_structure_results/PC3-4_bowhead_snps_LDprunedr08_n20_size1.png", width=6, height=5, dpi=1200)
#ggsave("rightwhale_map/pop_structure/PC3-4_bowhead_RWmap_snps_LDprunedr08_n20_size3.png", width=4, height=3, dpi=1200)
