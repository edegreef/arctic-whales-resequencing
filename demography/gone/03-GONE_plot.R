# Plot GONE results for a species at a time.

library(tidyverse)

# Let's do narwhal first
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/gone/with_maf05_10K_12chrscaf/")

# Load data
data1 <-  read.delim("Output_Ne_narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n57.12chrscaf.scafrename", header=T)

# Going to look at within last 200 generations
data1 <- subset(data1[1:200,])

# 21.9 for narwhal
gen=21.9

# Quick plot of the means (cut off at 200 generation time)
ggplot()+
  geom_line(data=data1, aes(x=(Generation),y=Geometric_mean), color="#2166AC", lwd=1.5)+
  theme_bw()+
  xlim(0,200)+
  xlab("generation")#+
  
  # and if want logscale:
  #scale_y_continuous(trans='log10')+
 # scale_x_continuous(trans='log10')+
 # annotation_logticks(sides = "lb")

## Using the 500 iterations to make a 95% confidence interval (based off of Kardos et al. 2023's script: https://github.com/martykardos/KillerWhaleInbreeding/blob/main/FigureCode/rCode_Fig_ED_1.R)

library(scales)
library(matrixStats)

# Load all the iteration files and put it in a matrix
files <- paste("outfileLD_TEMP/outfileLD_",1:500,"_GONE_Nebest",sep="")

# Create empty matrix
NeMat <- NULL

# Fill in matrix with all the info from the iterations
for(i in 1:500){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# Get CI for the recent 200 generations
NeCI <- matrix(NA,nrow=598,ncol=2)
for(i in 1:598){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

# Set up data to get ready for ggplot
NeCI_dat <- as.data.frame(NeCI)
NeCI_dat$Generation <- 1:nrow(NeCI_dat)

Ne_med <- as.data.frame(rowMedians(NeMat[1:500,]))
colnames(Ne_med) <- "median"
Ne_med$Generation <- 1:nrow(Ne_med)


# Plot
ggplot()+
  geom_ribbon(data=NeCI_dat, aes(x=Generation*gen, ymin=V1, ymax=V2), fill="#C6DBEF", alpha=0.5)+
  geom_line(data=Ne_med, aes(x=(Generation*gen),y=median), color="#6BAED6", lwd=1.5)+
  xlim(0,150*gen)+
  theme_bw()+
  xlab("Thousand years ago (kya)")+
  ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_y_continuous(breaks=c(0,10000,20000,30000,40000,50000),
                     labels=c("0e+00"=expression(0),"1e+04"=expression(10),"2e+04"=expression(20),"3e+04"=expression(30),"4e+04"=expression(40),"5e+04"=expression(50)))+
  scale_x_continuous(limits=c(0,150*gen),breaks=c(0,1000,2000,3000,4000),
                     labels=c("0"="0","1000"="1", "2000"="2","3000"="3","4000"="4"))

ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/plots/GONE_plot_narwhal_12chr_notlogged_150gen.png", height=3, width=3.5, dpi=1000) 


#########################################################################
# Bowhead whale time

setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/bowhead/gone/with_maf05/")

# Load data
data1 <-  read.delim("Output_Ne_bowhead_snps.filter1.miss.biallel.min100kb.autosomes.hwe.maf.n21.top200scaf", header=T)

# Going to look at within last 200 generations
data1 <- subset(data1[1:200,])

# 52.3 for bowhead
gen=52.3

# quick plot of the means (cut off at 200 generation time)
ggplot()+
  geom_line(data=data1, aes(x=(Generation),y=Geometric_mean), color="#2166AC", lwd=1.5)+
  theme_bw()+
  xlim(0,150)+
  xlab("generation")#+
  
  # and if want logscale:
#  scale_y_continuous(trans='log10')+
  #scale_x_continuous(trans='log10')+
#  annotation_logticks(sides = "lb")

## Using the 500 iterations to make a 95% confidence interval (based off of Kardos et al. 2023's script: https://github.com/martykardos/KillerWhaleInbreeding/blob/main/FigureCode/rCode_Fig_ED_1.R)

# Load all the iteration files and put it in a matrix
files <- paste("outfileLD_TEMP/outfileLD_",1:500,"_GONE_Nebest",sep="")

# Create empty matrix
NeMat <- NULL

# Fill in matrix with all the info from the iterations
for(i in 1:500){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# Get CI for the recent 200 generations
NeCI <- matrix(NA,nrow=598,ncol=2)
for(i in 1:598){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}


# Set up data to get ready for ggplot
NeCI_dat <- as.data.frame(NeCI)
NeCI_dat$Generation <- 1:nrow(NeCI_dat)

Ne_med <- as.data.frame(rowMedians(NeMat[1:500,]))
colnames(Ne_med) <- "median"
Ne_med$Generation <- 1:nrow(Ne_med)

# Plot
ggplot()+
  geom_ribbon(data=NeCI_dat, aes(x=Generation*gen, ymin=V1, ymax=V2), fill="#ffe38c", alpha=0.5)+
  geom_line(data=Ne_med, aes(x=(Generation*gen),y=median), color="#f8d568", lwd=1.5)+
  xlim(0,1500*gen)+
  theme_bw()+
  xlab("Thousand years ago (kya)")+
  ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_y_continuous(breaks=c(0,1e6,2e6,3e6,4e6,5e6),
                     labels=c("0e+00"=expression(0),"1e+06"="1,000","2e+06"="2,000","3e+06"="3,000","4e+06"="4,000","5e+06"="5,000"))+
  scale_x_continuous(limits=c(0,150*gen),breaks=c(0,2500,5000,7500,10000),
                     labels=c("0"="0","2500"="2.5", "5000"="5","7500"="7.5","10000"="10"))

ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/plots/GONE_plot_bowhead_top200_notlogged_150gen.png", height=3, width=3.5, dpi=1000) 

