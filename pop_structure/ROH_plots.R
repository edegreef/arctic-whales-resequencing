# Plotting ROH (runs of homozygosity)

# Modules
library(tidyverse)
library(patchwork)

# Starting with ROH
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/roh_redo_current_may2024/")

# Load ROH data for each pop
# Narwhal
NHB <- read.table("data/narwhal_NHB.ROH.nomaf.set17.hom.indiv", header=T)
EBB <- read.table("data/narwhal_EBB.ROH.nomaf.set17.hom.indiv", header=T)
WBB <- read.table("data/narwhal_WBB.ROH.nomaf.set17.hom.indiv", header=T)

# for bowhead (all as one pop)
BOW <- read.table("data/bowhead_RWmap.n20.ROH.nomaf.set17.hom.indiv", header=T)

# Mark region
BOW$region <- "BOW"
NHB$region <- "NHB"
EBB$region <- "EBB"
WBB$region <- "WBB"

# Merge files 
merge <- rbind(WBB, EBB, NHB, BOW)

set <- merge[, c("NSEG", "KB", "region")]

set$MB <- (set$KB)/1000

# Kind of by region- match the legend order
set$region <- factor(set$region , levels=c("WBB", "EBB", "NHB", "BOW")) 

#ggplot(set, aes(x=reorder(region,MB), y=MB, fill=region))+

bar_MB <- ggplot(set, aes(region, MB, fill=region))+
  #geom_boxplot(width=0.5, lwd=1,alpha=0.5)+
  geom_violin(width=1, lwd=0.8)+
  geom_boxplot(width=.15, fill="white", alpha=1)+
  theme_classic()+
  ylab("Total length of ROHs (Mb)")+
  xlab("Region")+
  scale_fill_manual(values=c("#F64F2E","#1EB4C4", "#807DBA", "orange"),  breaks=c("NHB", "EBB", "WBB", "BOW"))+
  theme(legend.position = "none")
bar_MB

#ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/roh_redo_current_may2024/ROH_violin.set8.png", width=4, height=4, dpi=600)


ggplot(data=set, aes(x=MB,y=NSEG))+
  geom_point(aes(colour=region), size=3, alpha=0.9)+
  theme_classic()+
  labs(color= "Region")+
  xlab("Total length of ROHs (Mb)")+
  ylab("Number of ROHs")+
  scale_colour_manual(name="Population",values = c("WBB"="#2171B5", "EBB"="#4292C6", "NHB"="#6BAED6", "BOW"="#f8d568"), labels=c("WBB"="Narwhal (WBB)", "EBB"="Narwhal (EBB)", "NHB"="Narwhal (NHB)", "BOW"="Bowhead whale"),breaks=c("BOW", "NHB", "EBB", "WBB"))
# scale_colour_manual(values=c("#E31A1C","#1D91C0", "#54278F"),breaks=c("RB", "MID", "NOR"))

#ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/roh_redo_current_may2024/ROH_points.set8.png", width=7, height=5, dpi=600)


set_bowhead_only <- subset(set, region=="BOW")
set_narwhal_only <- subset(set, region!="BOW")

scatter_narwhal <-ggplot(data=set_narwhal_only, aes(x=MB,y=NSEG))+
  annotation_custom(grid::linesGrob(gp = grid::gpar(col = "grey", lty = 2, lwd = 3))) +
  geom_point(aes(colour=region), size=3, alpha=0.9)+
  theme_classic()+
  labs(color= "Region")+
  xlab("Total length of ROHs (Mb)")+
  ylab("Number of ROHs")+
  scale_colour_manual(name="Population",values = c("WBB"="#2171B5", "EBB"="#6baed6", "NHB"="#bdd7e7", "BOW"="#f8d568"), labels=c("WBB"="Narwhal (WBB)", "EBB"="Narwhal (EBB)", "NHB"="Narwhal (NHB)", "BOW"="Bowhead whale"),breaks=c("BOW", "NHB", "EBB", "WBB"))+
  theme(legend.position = "none")

scatter_bowhead <-ggplot(data=set_bowhead_only, aes(x=MB,y=NSEG))+
  annotation_custom(grid::linesGrob(gp = grid::gpar(col = "grey", lty = 2, lwd = 3))) +
  geom_point(aes(colour=region), size=3, alpha=0.9)+
  theme_classic()+
  labs(color= "Region")+
  xlab("Total length of ROHs (Mb)")+
  ylab("Number of ROHs")+
  scale_colour_manual(name="Population",values = c("WBB"="#2171B5", "EBB"="#6baed6", "NHB"="#bdd7e7", "BOW"="#f8d568"), labels=c("WBB"="Narwhal (WBB)", "EBB"="Narwhal (EBB)", "NHB"="Narwhal (NHB)", "BOW"="Bowhead whale"),breaks=c("BOW", "NHB", "EBB", "WBB"))+
  theme(legend.position = "none")


scatter_narwhal + scatter_bowhead

#ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/roh_redo_current_may2024/ROH_scatter_sepsp_set8_maf_RWmap.png", height=3.5, width=10, dpi=900)

scatter_narwhal
#ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/roh_redo_current_may2024/ROH_scatter_narwhal_set8_NOmaf_DPI.png", height=3, width=3, dpi=1200)
scatter_bowhead
#ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/roh_redo_current_may2024/ROH_scatter_bowhead_RWmap_set8_NOmaf_DPI.png", height=3, width=3, dpi=1200)


######################

### Next going to look into FROH

# Adding extra column to data for FROH - fraction of genome as ROH
scafs_bow <- read.table("C:/Users/eveli/Dropbox/Whales with Garroway/ref_genomes/right_whale/GCA_028564815.2_mEubGla1.1.hap2._XY_genomic.fna.fai", header=F)
colnames(scafs_bow) <- c("scaffold", "length", "V3", "V4", "V5")
scafs_nar <- read.table("C:/Users/eveli/Dropbox/Whales with Garroway/ref_genomes/narwhal/NAR_GCF_005190385.2.scafname.fasta.fai", header=F)
colnames(scafs_nar) <- c("scaffold", "length", "V3", "V4", "V5")


# remember to just the scaffolds > 100kb (and only autosomes)
scafs_bow_100kb <- subset(scafs_bow, length >= 100000)
scafs_nar_100kb <- subset(scafs_nar, length >= 100000)

# sex chr bowhead
#--not-chr NC_045806.1 --not-chr NC_045807.1
scafs_bow_100kb_auto <- subset(scafs_bow_100kb, scaffold != "CM053059.1" & scaffold != "CM053080.1")

total_length_bow_100kb_auto <- sum(scafs_bow_100kb_auto$length)
total_length_bow_100kb_auto
# 2677992299

# sex chr narwhal
x_chr_n <- read.table("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/narwhal/difcover/final_scaffolds_Xlinked.txt", header=F)
colnames(x_chr_n) <- "scaffold"
y_chr_n <- read.table("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/narwhal/difcover/final_scaffolds_Ylinked.txt", header=F)
colnames(y_chr_n) <- "scaffold"

scafs_nar_100kb_auto <- anti_join(scafs_nar_100kb, x_chr_n, by="scaffold")
scafs_nar_100kb_auto <- anti_join(scafs_nar_100kb, y_chr_n, by="scaffold")

total_length_nar_100kb_auto <- sum(scafs_nar_100kb_auto$length)
total_length_nar_100kb_auto
# 2322249225


# Remember to make it Mb unit to be consistent
# 1 million bp per Mb
total_length_bow_100kb_auto_MB <- total_length_bow_100kb_auto/1000000
total_length_bow_100kb_auto_MB
#2677.992

total_length_nar_100kb_auto_MB <- total_length_nar_100kb_auto/1000000
total_length_nar_100kb_auto_MB
#2322.249


# Estimate FROH
# redo the merge to make sure sp with respective genome
merge_nar <- rbind(WBB, EBB, NHB)
merge_bow <- BOW

set_nar <- merge_nar[, c("NSEG", "KB", "region")]
set_bow <- merge_bow[, c("NSEG", "KB", "region")]

set_nar$MB <- (set_nar$KB)/1000
set_bow$MB <- (set_bow$KB)/1000

set_nar$Fall <- (set_nar$MB/total_length_nar_100kb_auto_MB)
set_bow$Fall <- (set_bow$MB/total_length_bow_100kb_auto_MB)

#now merge
set_final <- rbind(set_nar, set_bow)

set_final$row <- 1:nrow(set_final)

# Plot it- going to stick with the autosomes only b/c that's how data was mapped
Froh <- ggplot(data=set_final, aes(x=reorder(row,Fall),y=Fall))+
  theme_bw()+
  geom_hline(yintercept=0.01, color="grey88")+
  geom_hline(yintercept=0.02, color="grey88")+
  geom_hline(yintercept=0.03, color="grey88")+
  geom_hline(yintercept=0.04, color="grey88")+
  geom_hline(yintercept=0.05, color="grey88")+
  
 # geom_hline(yintercept=0.06, color="grey88")+
  labs(fill= "Pop")+
  xlab("Individual")+
  ylab(expression(paste(italic("F")["ROH"],sep="")))+
  # scale_fill_manual(values = c("P2"="#D6604D","P1"="#4393C3"), labels=(c("P2"="P2", "P1"="P1")))+
  scale_fill_manual(name="Population",values = c("WBB"="#2171B5", "EBB"="#4292C6", "NHB"="#6BAED6", "BOW"="#f8d568"), labels=c("WBB"="Narwhal (WBB)", "EBB"="Narwhal (EBB)", "NHB"="Narwhal (NHB)", "BOW"="Bowhead whale"),breaks=c("BOW", "NHB", "EBB", "WBB"))+
  
  geom_bar(stat = "identity", aes(fill=region), colour="black", linewidth=0.15, width=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x=element_blank())+
  scale_y_continuous(expand=c(0,0))#, limits=c(0,0.07))

Froh
#ggsave("froh_plot_maf_RWmap_set8.png", width=7, height=2.5, dpi=800)


######################################
# To use ROH > ## Mb then look at the .hom file instead of .hom.indiv.
##### UPDATE >>> 1.5Mb is too big for bowhead, so going to use 100kb
# load files
# Load ROH data for each pop

NHB <- read.table("data/narwhal_NHB.ROH.nomaf.set19.hom", header=T)
EBB <- read.table("data/narwhal_EBB.ROH.nomaf.set19.hom", header=T)
WBB <- read.table("data/narwhal_WBB.ROH.nomaf.set19.hom", header=T)

# for bowhead (all as one pop)
BOW <- read.table("data/bowhead_RWmap.n20.ROH.nomaf.set19.hom", header=T)
#BOW <- subset(BOW, FID!="RMD_BM_96_1")
#BOW <- subset(BOW, FID!="AR_BM_SH_2003_01")

# Filter at ## KB min length
WBB_hom <- subset(WBB, KB > 100)
EBB_hom <- subset(EBB, KB > 100)
NHB_hom <- subset(NHB, KB > 100)
BOW_hom <- subset(BOW, KB > 100)

# Then pull up sums by individual (like making own .hom.indiv file)
WBB_hom_indiv <- WBB_hom %>% group_by(FID) %>% summarize(total_kb = sum(KB))
EBB_hom_indiv <- EBB_hom %>% group_by(FID) %>% summarize(total_kb = sum(KB))
NHB_hom_indiv <- NHB_hom %>% group_by(FID) %>% summarize(total_kb = sum(KB))
BOW_hom_indiv <- BOW_hom %>% group_by(FID) %>% summarize(total_kb = sum(KB))

# Use total_length_10Mb_auto_MB as length of autosomes in this dataset.
# Put everything in MB to keep consistent
# Total_length_10Mb_auto_MB
WBB_hom_indiv$total_MB <- WBB_hom_indiv$total_kb/1000
EBB_hom_indiv$total_MB <- EBB_hom_indiv$total_kb/1000
NHB_hom_indiv$total_MB <- NHB_hom_indiv$total_kb/1000
BOW_hom_indiv$total_MB <- BOW_hom_indiv$total_kb/1000

WBB_hom_indiv$pop <- "WBB"
EBB_hom_indiv$pop <- "EBB"
NHB_hom_indiv$pop <- "NHB"
BOW_hom_indiv$pop <- "BOW"

# Combine into one
set_nar <- rbind(WBB_hom_indiv, EBB_hom_indiv, NHB_hom_indiv)
set_bow <- BOW_hom_indiv

# Add the proportion (FROH)
set_nar$FROH <- (set_nar$total_MB/total_length_nar_100kb_auto_MB)
set_bow$FROH <- (set_bow$total_MB/total_length_bow_100kb_auto_MB)

set_final <- rbind(set_nar, set_bow)

set_final$row <- 1:nrow(set_final)

# Plot it
Froh <- ggplot(data=set_final, aes(x=reorder(row,FROH),y=FROH))+
  theme_bw()+
  geom_hline(yintercept=0.04, color="grey88")+
 # geom_hline(yintercept=0.06, color="grey88")+
  labs(fill= "Pop")+
  xlab("Individual")+
  ylab(expression(paste(italic("F")["ROH(>100kb)"],sep="")))+
  scale_fill_manual(name="Population",values = c("WBB"="#2171B5", "EBB"="#4292C6", "NHB"="#6BAED6", "BOW"="#f8d568"), labels=c("WBB"="Narwhal (WBB)", "EBB"="Narwhal (EBB)", "NHB"="Narwhal (NHB)", "BOW"="Bowhead whale"),breaks=c("BOW", "NHB", "EBB", "WBB"))+
  geom_bar(stat = "identity", aes(fill=pop), colour="black", linewidth=0.15, width=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x=element_blank())+
  scale_y_continuous(expand=c(0,0))#, limits=c(0,0.07))

Froh

# Plot as box plot instead
ggplot(data=set_final, aes(x=pop,y=FROH))+
  theme_bw()+
 # geom_hline(yintercept=0.04, color="grey88")+
  # geom_hline(yintercept=0.06, color="grey88")+
  labs(fill= "Pop")+
  xlab("Individual")+
  ylab(expression(paste(italic("F")["ROH"],sep="")))+
  scale_color_manual(name="Population",values = c("WBB"="#2171B5", "EBB"="#6baed6", "NHB"="#bdd7e7", "BOW"="#f8d568"), labels=c("WBB"="Narwhal (WBB)", "EBB"="Narwhal (EBB)", "NHB"="Narwhal (NHB)", "BOW"="Bowhead whale"),breaks=c("BOW", "NHB", "EBB", "WBB"))+
 # geom_violin(width=0.5, lwd=1, fill=NA, colour="gray70")+
  #geom_violin(width=1, lwd=0.8, aes(fill=pop))+
 #  geom_jitter(width=0.15, aes(fill=pop), alpha=0.7)+
#geom_boxplot(width=.15, fill="white", alpha=1)+
  geom_boxplot(width=0.5, lwd=1, fill=NA, colour="gray65")+
  geom_jitter(width=0.15, aes(color=pop),size=2, alpha=0.9)+
  
 # geom_bar(stat = "identity", aes(fill=pop), colour="black", linewidth=0.15, width=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x=element_blank(),
        axis.title=element_text(size=14))+
  scale_x_discrete(limits=c("NHB", "EBB", "WBB", "BOW"))+
  scale_y_continuous(limits=c(0,0.05))


#ggsave("froh_boxplot_maf_set8_RWmap_size2.png", width=5.5, height=5, dpi=800)
#ggsave("froh_boxplot_maf_set8_RWmap_size3.png", width=8, height=5, dpi=800)

