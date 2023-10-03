# Plotting Tajd, pi, and ROH (runs of homozygosity)

# Modules
library(ggplot2)
library(patchwork)
#library(ggbeeswarm)

# Starting with ROH
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/narwhal/taj_pi_roh")

# for bowhead (all as one pop)
#BOW <- read.table("bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.n21.hom.indiv", header=T)

#(AB, BI, CR, GF, IG, PB, PG, PI, RB, RE, SB)

# Load ROH data for each site
AB <- read.table("AB.hom.indiv", header=T)
BI <- read.table("BI.hom.indiv", header=T)
CR <- read.table("CR.hom.indiv", header=T)
GF <- read.table("GF.hom.indiv", header=T)
IG <- read.table("IG.hom.indiv", header=T)
PB <- read.table("PB.hom.indiv", header=T)
#PG <- read.table("PG.hom.indiv", header=T)
PI <- read.table("PI.hom.indiv", header=T)
RB <- read.table("RB.hom.indiv", header=T)
RE <- read.table("RE.hom.indiv", header=T)
SB <- read.table("SB.hom.indiv", header=T)

# And for subgroups
MID <- read.table("MID.hom.indiv", header=T)
NOR <- read.table("NOR.hom.indiv", header=T)

# Bowhead as one
BOW$region <- "BOW"
sample_info <- read.csv("C:/Users/eveli/Dropbox/Whales with Garroway/bowhead/sample_info_bowhead.csv")
# for bowhead:
sample_info <- subset(sample_info, sample_ID!="99_01") #kin pair B. for n22
sample_info <- subset(sample_info, sample_ID!="ARBMGH_2002_001") #kin pair A. for n21

BOW$region <- sample_info$location_ID
BOW$region <- "BOW"

# Narwhal
AB$region <- "AB"
BI$region <- "BI"
CR$region <- "CR"
GF$region <- "GF"
IG$region <- "IG"
PB$region <- "PB"
#PG$region <- "PG"
PI$region <- "PI"
RB$region <- "NHB"
RE$region <- "RE"
SB$region <- "SB"

MID$region <- "EBB"
NOR$region <- "WBB"

# Merge files 
#merge <- BOW
merge <- rbind(AB, BI, CR, GF, PB, PI, RB, RE, SB, IG) #exclude PG
merge <- rbind(NOR, MID, RB) #exclude PG

set <- merge[, c("NSEG", "KB", "region")]

set$MB <- (set$KB)/1000

roh <- ggplot(data=set, aes(x=NSEG,y=MB))+
  geom_point(aes(colour=region), size=3, alpha=0.9)+
  theme_classic()+
  labs(color= "Region")+
  xlab("Total length of ROHs (Mb)")+
  ylab("Number of ROHs")+
  scale_colour_manual(values=c("#54278F", "#807DBA", "#BCBDDC", "#0C2C84", "#225EA8", "#1D91C0", "#41B6C4", "#7FCDBB", "#FED976", "#FD8D3C", "#E31A1C"),             
                      breaks=c("GF", "RE", "SB", "IG", "PB", "AB", "PI", "CR", "BI", "PG", "RB"))
 # scale_colour_manual(values=c("#E31A1C","#1D91C0", "#54278F"),breaks=c("RB", "MID", "NOR"))
roh

#ggsave("roh_linear_plot_colour_3group.png", width=5, height=4, dpi=700)
# trying to change order


# Kind of by region- match the legend order
set$region <- factor(set$region , levels=c("GF", "RE", "SB", "PB", "IG", "AB", "PI", "CR", "BI", "RB")) 

set$region <- factor(set$region , levels=c("WBB", "EBB", "NHB")) 

#ggplot(set, aes(x=reorder(region,MB), y=MB, fill=region))+

bar_MB <- ggplot(set, aes(region, MB, fill=region))+
  #geom_boxplot(width=0.5, lwd=1,alpha=0.5)+
 geom_violin(width=1, lwd=0.8)+
  geom_boxplot(width=.15, fill="white", alpha=1)+
  
   # geom_boxplot(width=0.7, lwd=1, fill=NA, colour="gray88")+
 # geom_jitter(width=0.25,pch=21, colour="black",  size=4)+
  #geom_point(shape = 21, colour = "black", fill = "gray",size=1)+
 # geom_quasirandom(pch=21,size=2, colour="black", width=0.2, show.legend = FALSE, alpha =1)+
  theme_classic()+
  ylab("Total length of ROHs (Mb)")+
  xlab("Region")+
 # ggtitle("ROH")+
 # scale_fill_manual(values=c("#54278F", "#807DBA", "#BCBDDC", "#0C2C84", "#225EA8", "#1D91C0", "#41B6C4", "#7FCDBB", "#FED976", "#FD8D3C", "#E31A1C"),             
                #      breaks=c("GF", "RE", "SB", "IG", "PB", "AB", "PI", "CR", "BI", "PG", "RB"))+
  scale_fill_manual(values=c("#F64F2E","#1EB4C4", "#807DBA"),  breaks=c("NHB", "EBB", "WBB"))+
  theme(legend.position = "none")
bar_MB
#ggsave("ROH_3groups_violin_BBupdated.png", width=3.5, height=3, dpi=600)

# TAJ
AB <- read.table("AB.10kb.Tajima.D", header=T)
BI <- read.table("BI.10kb.Tajima.D", header=T)
CR <- read.table("CR.10kb.Tajima.D", header=T)
GF <- read.table("GF.10kb.Tajima.D", header=T)
IG <- read.table("IG.10kb.Tajima.D", header=T)
PB <- read.table("PB.10kb.Tajima.D", header=T)
#PG <- read.table("PG.5kb.Tajima.D", header=T)
PI <- read.table("PI.10kb.Tajima.D", header=T)
RB <- read.table("RB.10kb.Tajima.D", header=T)
RE <- read.table("RE.10kb.Tajima.D", header=T)
SB <- read.table("SB.10kb.Tajima.D", header=T)

MID <- read.table("MID.10kb.Tajima.D", header=T)
NOR <- read.table("NOR.10kb.Tajima.D", header=T)

AB$region <- "AB"
BI$region <- "BI"
CR$region <- "CR"
GF$region <- "GF"
IG$region <- "IG"
PB$region <- "PB"
#PG$region <- "PG"
PI$region <- "PI"
RB$region <- "RB"
RE$region <- "RE"
SB$region <- "SB"

MID$region <- "MID"
NOR$region <- "NOR"

merge_taj <- rbind(AB, BI, CR, GF,PB, PI, RB, RE, SB, IG) # exclude PG
merge_taj <- rbind(NOR, MID, RB)

# kind of by region- match the legend order?
merge_taj$region <- factor(merge_taj$region , levels=c("GF", "RE", "SB", "PB", "IG", "AB", "PI", "CR", "BI", "RB")) 
# kind of by region- match the legend order?
merge_taj$region <- factor(merge_taj$region , levels=c("NOR", "MID", "RB"))
                           
taj <- ggplot(merge_taj, aes(region, TajimaD, fill=region))+
 # ggplot(merge_taj, aes(x=reorder(region,TajimaD,na.rm=TRUE), y=TajimaD, fill=region))+
  #geom_boxplot(width=0.5, lwd=1, fill=NA, color="gray88")+
  
  #geom_boxplot(width=0.5, lwd=0.8)+
  geom_violin(width=0.7, lwd=0.8)+
  geom_boxplot(width=.15, fill="white", alpha=1)+
  
  #geom_boxplot(width=0.7, lwd=1, fill=NA, colour="gray88")+
 # geom_jitter(width=0.05,pch=21, colour="black",  size=2, alpha=0.5)+
  #geom_quasirandom(pch=21,size=2, colour="black", width=0.2, show.legend = FALSE, alpha =0.5)+
  theme_classic()+
  ylab("Tajima's D")+
  xlab("Region")+
 # ggtitle("Tajima's D")+
  theme(legend.position="none")+
  #scale_fill_manual(values=c("#54278F", "#807DBA", "#BCBDDC", "#0C2C84", "#225EA8", "#1D91C0", "#41B6C4", "#7FCDBB", "#FED976", "#FD8D3C", "#E31A1C"),             
                    #  breaks=c("GF", "RE", "SB", "IG", "PB", "AB", "PI", "CR", "BI", "PG", "RB"))+
  scale_fill_manual(values=c("#F64F2E","#1EB4C4", "#807DBA"),  breaks=c("RB", "MID", "NOR"))+
  geom_hline(yintercept=0, colour="gray")
taj


# PI
AB <- read.table("AB.10kb.windowed.pi", header=T)
BI <- read.table("BI.10kb.windowed.pi", header=T)
CR <- read.table("CR.10kb.windowed.pi", header=T)
GF <- read.table("GF.10kb.windowed.pi", header=T)
IG <- read.table("IG.10kb.windowed.pi", header=T)
PB <- read.table("PB.10kb.windowed.pi", header=T)
#PG <- read.table("PG.5kb.windowed.pi", header=T)
PI <- read.table("PI.10kb.windowed.pi", header=T)
RB <- read.table("RB.10kb.windowed.pi", header=T)
RE <- read.table("RE.10kb.windowed.pi", header=T)
SB <- read.table("SB.10kb.windowed.pi", header=T)

MID <- read.table("MID.10kb.windowed.pi", header=T)
NOR <- read.table("NOR.10kb.windowed.pi", header=T)


AB$region <- "AB"
BI$region <- "BI"
CR$region <- "CR"
GF$region <- "GF"
IG$region <- "IG"
PB$region <- "PB"
#PG$region <- "PG"
PI$region <- "PI"
RB$region <- "RB"
RE$region <- "RE"
SB$region <- "SB"

MID$region <- "MID"
NOR$region <- "NOR"


merge_pi <- rbind(AB, BI, CR, GF, PB, PI, RB, RE, SB, IG) #exclude PG
merge_pi <- rbind(NOR, MID, RB) #exclude PG


# kind of by region- match the legend order?
merge_pi$region <- factor(merge_pi$region , levels=c("GF", "RE", "SB", "PB", "IG", "AB", "PI", "CR", "BI", "RB")) 
merge_pi$region <- factor(merge_pi$region , levels=c("NOR", "MID", "RB"))
  
pi <- ggplot(merge_pi, aes(region, PI, fill=region))+
  #ggplot(merge_pi, aes(x=reorder(region,PI,na.rm=TRUE), y=PI, fill=region))+
 # geom_boxplot(width=0.5, lwd=0.8)+
  geom_violin(width=0.7, lwd=0.8)+
  geom_boxplot(width=.15, fill="white", alpha=1)+
 # geom_boxplot(width=0.7, lwd=1, fill=NA, colour="gray88")+
  # geom_jitter(width=0.05, colour="gray50", size=2)+
#  geom_jitter(width=0.05,pch=21, colour="black",  size=2, alpha=0.5)+
  theme_classic()+
  ylab("Pi")+
  xlab("Region")+
#  ggtitle("Nucleotide Diversity")+
 # scale_fill_manual(values=c("#54278F", "#807DBA", "#BCBDDC", "#0C2C84", "#225EA8", "#1D91C0", "#41B6C4", "#7FCDBB", "#FED976", "#FD8D3C", "#E31A1C"),             
    #                breaks=c("GF", "RE", "SB", "IG", "PB", "AB", "PI", "CR", "BI", "PG", "RB"))
scale_fill_manual(values=c("#F64F2E","#1EB4C4", "#807DBA"),  breaks=c("RB", "MID", "NOR"))+
  theme(legend.position="none")
pi

bar_MB / taj / pi

taj+pi

ggsave("taj10kb_pi10kb_roh_plot_omitPG.png", width=6, height=8, dpi=1000)
ggsave("taj10kb_pi10kb_plot_3group.png", width=8, height=4, dpi=1000)
