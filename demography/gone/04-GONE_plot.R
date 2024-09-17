# plotting GONE results for a species at a time.

library(tidyverse)

## Using the 500 iterations to make a 95% confidence interval (based off of Kardos et al. 2023's script: https://github.com/martykardos/KillerWhaleInbreeding/blob/main/FigureCode/rCode_Fig_ED_1.R)

library(scales)
library(matrixStats)


####### Bowhead
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/bowhead/rightwhale_map/gone/param_original/bowhead_RWmap_10MB_nomaf_500iter")

data <-  read.delim("Output_Ne_bowhead_RWmap_snps_rmfix.filter1.miss.biallel.min100kb.autosomes.hwe.n20.10MB.scafrename", header=T)

# going to look at within last 150 generations
data <- subset(data[1:150,])

# quick plot of the means
ggplot()+
  geom_line(data=data, aes(x=(Generation),y=Geometric_mean), color="#2166AC", lwd=1.5)+
  theme_bw()+
  xlim(0,150)+
  xlab("generation")+
  
  # and if want logscale:
 scale_y_continuous(trans='log10')+
 scale_x_continuous(trans='log10')+
  annotation_logticks(sides = "lb")

# load all the iteration files and put it in a matrix
files <- paste("outfileLD_TEMP/outfileLD_",1:500,"_GONE_Nebest",sep="")
# create empty matrix
NeMat <- NULL

# fill in matrix with all the info from the iterations
for(i in 1:500){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# get CI for the recent 150 generations
NeCI <- matrix(NA,nrow=789,ncol=2)
for(i in 1:789){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

# set up data to get ready for ggplot
NeCI_dat <- as.data.frame(NeCI)
NeCI_dat$Generation <- 1:nrow(NeCI_dat)

Ne_med <- as.data.frame(rowMedians(NeMat[1:789,]))
colnames(Ne_med) <- "median"
Ne_med$Generation <- 1:nrow(Ne_med)

BOW <- ggplot()+
  geom_ribbon(data=NeCI_dat, aes(x=Generation, ymin=V1, ymax=V2), fill="#ffe38c", alpha=0.5)+
  geom_line(data=Ne_med, aes(x=(Generation),y=median), color="#f8d568", lwd=1.5)+
  xlim(0,150)+
  theme_bw()+
  xlab("Generations")+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  #ylab(expression(paste("Effective population size")))
  ylab(expression(paste("Effective population size (")*10^{3}*")"))+
  scale_y_continuous(trans='log10', breaks=c(1000,10000,100000), limits=c(500,200000), labels=c("1e+03"=expression(1),"1e+04"=expression(10), "1e+05"=expression(100)), expand=c(0,0))+
  scale_x_continuous(trans='log10', limits =c(1,150), expand=c(0,0))+
  annotation_logticks(sides = "lb")

BOW

ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/plots/GONE-revised-bowhead_RWmap_10MB_nomaf_500iter_originalparam.png", width=4, height=3, dpi=800)

####### Narwhal
#colors
"WBB"="#2171B5", "EBB"="#6baed6", "NHB"="#bdd7e7"

## EBB (EBB was old abbreviation, noting here this is for Baffin Island)
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/gone/param_original/nomaf_10MB_EBB_500iter/")

data_EBB <-  read.delim("Output_Ne_narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.EBB.10MB.scafrename", header=T)

# going to look at within last 150 generations
data_EBB <- subset(data_EBB[1:150,])

# quick plot of the means
ggplot()+
  geom_line(data=data_EBB, aes(x=(Generation),y=Geometric_mean), color="#2166AC", lwd=1.5)+
  theme_bw()+
  xlim(0,150)+
  xlab("generation")

# and if want logscale:
#scale_y_continuous(trans='log10')+
#scale_x_continuous(trans='log10')+
#annotation_logticks(sides = "lb")

# load all the iteration files and put it in a matrix
files <- paste("outfileLD_TEMP/outfileLD_",1:500,"_GONE_Nebest",sep="")

# create empty matrix
NeMat <- NULL

# fill in matrix with all the info from the iterations
for(i in 1:500){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# get CI for the recent 200 generations
NeCI <- matrix(NA,nrow=789,ncol=2)
for(i in 1:789){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}


# set up data to get ready for ggplot
NeCI_dat_EBB <- as.data.frame(NeCI)
NeCI_dat_EBB$Generation <- 1:nrow(NeCI_dat_EBB)

Ne_med_EBB <- as.data.frame(rowMedians(NeMat[1:789,]))
colnames(Ne_med_EBB) <- "median"
Ne_med_EBB$Generation <- 1:nrow(Ne_med_EBB)

EBB_plot <- ggplot()+
  geom_ribbon(data=NeCI_dat_EBB, aes(x=Generation, ymin=V1, ymax=V2), fill="#6BAED6", alpha=0.5)+
  geom_line(data=Ne_med_EBB, aes(x=(Generation),y=median), color="#2171B5", lwd=1.5)+
  xlim(0,150)+
  theme_bw()+
  xlab("Generations")+
 # ylab("Effective population size")+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  ylab(expression(paste("Effective population size (")*10^{3}*")"))+
  scale_y_continuous(trans='log10', breaks=c(10000,100000), limits=c(1000,1000000), labels=c("1e+04"=expression(10), "1e+05"=expression(100)), expand=c(0,0))+
  scale_x_continuous(trans='log10', limits =c(1,150), expand=c(0,0))+
  annotation_logticks(sides = "lb")
EBB_plot


#ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/plots/GONE-revised-narwhal_EBB_10MB_nomaf_500iter_originalparam.png", width=8, height=3, dpi=800)

## WBB (WBB was old abbreviation, noting here this is for Canadian High Arctic)
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/gone/param_original/nomaf_10MB_WBB_500iter/")

data_WBB <-  read.delim("Output_Ne_narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.WBB.10MB.scafrename", header=T)

# going to look at within last 150 generations
data_WBB <- subset(data_WBB[1:150,])

# quick plot of the means
ggplot()+
  geom_line(data=data_EBB, aes(x=(Generation),y=Geometric_mean), color="#2166AC", lwd=1.5)+
  theme_bw()+
  xlim(0,150)+
  xlab("generation")

# and if want logscale:
#scale_y_continuous(trans='log10')+
#scale_x_continuous(trans='log10')+
#annotation_logticks(sides = "lb")

# load all the iteration files and put it in a matrix
files <- paste("outfileLD_TEMP/outfileLD_",1:500,"_GONE_Nebest",sep="")

# create empty matrix
NeMat <- NULL

# fill in matrix with all the info from the iterations
for(i in 1:500){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# get CI for the recent 200 generations
NeCI <- matrix(NA,nrow=789,ncol=2)
for(i in 1:789){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}


# set up data to get ready for ggplot
NeCI_dat_WBB <- as.data.frame(NeCI)
NeCI_dat_WBB$Generation <- 1:nrow(NeCI_dat_WBB)

Ne_med_WBB <- as.data.frame(rowMedians(NeMat[1:789,]))
colnames(Ne_med_WBB) <- "median"
Ne_med_WBB$Generation <- 1:nrow(Ne_med_WBB)

WBB_plot <- ggplot()+
  geom_ribbon(data=NeCI_dat_WBB, aes(x=Generation, ymin=V1, ymax=V2), fill="#2171B5", alpha=0.5)+
  geom_line(data=Ne_med_WBB, aes(x=(Generation),y=median), color="#08306B", lwd=1.5)+
  xlim(0,150)+
  theme_bw()+
  xlab("Generations")+
  # ylab("Effective population size")+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  ylab(expression(paste("Effective population size (")*10^{3}*")"))+
  scale_y_continuous(trans='log10', breaks=c(10000,100000), limits=c(1000,1000000), labels=c("1e+04"=expression(10), "1e+05"=expression(100)), expand=c(0,0))+
  scale_x_continuous(trans='log10', limits =c(1,150), expand=c(0,0))+
  annotation_logticks(sides = "lb")
WBB_plot

## RB (RB was old abbreviation, noting here this is for Northern Hudson Bay)
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/gone/param_original/nomaf_10MB_RB_500iter/")

data_RB <-  read.delim("Output_Ne_narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.RB.10MB.scafrename", header=T)

# going to look at within last 150 generations
data_RB <- subset(data_RB[1:150,])

# quick plot of the means
ggplot()+
  geom_line(data=data_RB, aes(x=(Generation),y=Geometric_mean), color="#2166AC", lwd=1.5)+
  theme_bw()+
  xlim(0,150)+
  xlab("generation")+
  
  # and if want logscale:
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  annotation_logticks(sides = "lb")

# load all the iteration files and put it in a matrix
files <- paste("outfileLD_TEMP/outfileLD_",1:500,"_GONE_Nebest",sep="")

# create empty matrix
NeMat <- NULL

# fill in matrix with all the info from the iterations
for(i in 1:500){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# get CI for the recent 200 generations
NeCI <- matrix(NA,nrow=789,ncol=2)
for(i in 1:789){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}


# set up data to get ready for ggplot
NeCI_dat_RB <- as.data.frame(NeCI)
NeCI_dat_RB$Generation <- 1:nrow(NeCI_dat_RB)

Ne_med_RB <- as.data.frame(rowMedians(NeMat[1:789,]))
colnames(Ne_med_RB) <- "median"
Ne_med_RB$Generation <- 1:nrow(Ne_med_RB)

RB_plot <- ggplot()+
  geom_ribbon(data=NeCI_dat_RB, aes(x=Generation, ymin=V1, ymax=V2), fill="#DEEBF7", alpha=0.5)+
  geom_line(data=Ne_med_RB, aes(x=(Generation),y=median), color="#9ECAE1", lwd=1.5)+
  xlim(0,150)+
  theme_bw()+
  xlab("Generations")+
  ylab(expression(paste("Effective population size (")*10^{3}*")"))+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_y_continuous(trans='log10', breaks=c(10000,100000), limits=c(1000,1000000), labels=c("1e+04"=expression(10), "1e+05"=expression(100)), expand=c(0,0))+
  scale_x_continuous(trans='log10', limits =c(1,150), expand=c(0,0))+
  annotation_logticks(sides = "lb")
RB_plot


library(patchwork)
(WBB_plot + EBB_plot) / (RB_plot + BOW)

# Plot with narwhal subgroups together
narwhal_all <- ggplot()+
  geom_ribbon(data=NeCI_dat_EBB, aes(x=Generation, ymin=V1, ymax=V2), fill="#6BAED6", alpha=0.5)+
  geom_line(data=Ne_med_EBB, aes(x=(Generation),y=median), color="#2171B5", lwd=1.5)+
  geom_ribbon(data=NeCI_dat_WBB, aes(x=Generation, ymin=V1, ymax=V2), fill="#2171B5", alpha=0.5)+
  geom_line(data=Ne_med_WBB, aes(x=(Generation),y=median), color="#08306B", lwd=1.5)+
  geom_ribbon(data=NeCI_dat_RB, aes(x=Generation, ymin=V1, ymax=V2), fill="#DEEBF7", alpha=0.5)+
  geom_line(data=Ne_med_RB, aes(x=(Generation),y=median), color="#9ECAE1", lwd=1.5)+
  xlim(0,150)+
  theme_bw()+
  xlab("Generations")+  
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
 #ylab("Effective population size")
  ylab(expression(paste("Effective population size (")*10^{3}*")"))+
  scale_y_continuous(trans='log10', breaks=c(10000,100000), limits=c(1000,1000000), labels=c("1e+04"=expression(10), "1e+05"=expression(100)), expand=c(0,0))+
  scale_x_continuous(trans='log10', limits =c(1,150), expand=c(0,0))+
  annotation_logticks(sides = "lb")
narwhal_all

ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/plots/GONE-revised-narwhal_subgroup_together_10MB_nomaf_500iter_originalparam.png", width=4, height=3, dpi=800)
