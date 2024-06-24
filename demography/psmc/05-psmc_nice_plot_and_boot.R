# Plot PSMC data as main run and also with bootstrapping
# First part for Narwhal, second part for Bowhead whale

library(tidyverse)

############# NARWHAL
# Quick plot without bootstrap
# EBB
psmc_dat_EBB <- read.delim("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/psmc/ARCR_07_1065/ARCR_07_1065_1.56e-8_g21.9plot.0.txt", header=FALSE)
#WBB
psmc_dat_WBB <- read.delim("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/psmc/ARGF_07_1127/ARGF_07_1127_1.56e-8_g21.9plot.0.txt", header=FALSE)
#NHB
psmc_dat_NHB <- read.delim("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/psmc/ARRB_99_1026/ARRB_99_1026_1.56e-8_g21.9plot.0.txt", header=FALSE)
  
# plotting just 1
ggplot() + 
  geom_step(data=psmc_dat_EBB, aes(x=V1, y=V2), color="#1EB4C4") + 
  geom_step(data=psmc_dat_WBB, aes(x=V1, y=V2), color="#807DBA") + 
  geom_step(data=psmc_dat_NHB, aes(x=V1, y=V2), color="#F64F2E") + 
  scale_x_continuous(trans='log10') + 
  annotation_logticks(sides = "lb") +
  theme_classic()


# Then load in bootstrap files (change dir first)
# create list of file names to load (100 per group is too much to load individually)
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/psmc/ARCR_07_1065/boot_mu1.56")
boots_filenames_EBB <- list.files(pattern="*.txt")

# read data in
boots_EBB <- lapply(boots_filenames_EBB,function(i){
  read.delim(i, header=FALSE)
})

boots_all_EBB <- do.call(rbind, boots_EBB)

# every 58 rows is new boot. need column to label this as run number.
boots_all_EBB$run <- rep(c(1:101), each=58)
boots_all_EBB$pop <- "EBB"

# merge
runs_EBB <- rbind(boots_all_EBB)

# add unique id for label+run
runs_EBB$plot_id <- paste(runs_EBB$pop, runs_EBB$run, sep="_")

# adjusting order to make legend by latitude
runs_EBB$pop <- factor(runs_EBB$pop, levels = c("EBB"))

# extract median run per group
psmc_boot_EBB <- boots_all_EBB %>% group_by(run) %>% mutate(row=row_number())
psmc_boot_med_EBB <- aggregate(psmc_boot_EBB[,c("V1","V2")], list(psmc_boot_EBB$row), FUN=median)
colnames(psmc_boot_med_EBB)[1] <- "row"
psmc_boot_med_EBB$pop <- "EBB"


# merge
MED_EBB <- rbind(psmc_boot_med_EBB)

# adjust label order by latitude (for legend)
MED_EBB$pop <- factor(MED_EBB$pop, levels = c("EBB"))

########### same thing again but for WBB - should probably make this into a function for efficiency ...
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/psmc/ARGF_07_1127/boot_mu1.56")
boots_filenames_WBB <- list.files(pattern="*.txt")
boots_WBB <- lapply(boots_filenames_WBB,function(i){
  read.delim(i, header=FALSE)
})
boots_all_WBB <- do.call(rbind, boots_WBB)
boots_all_WBB$run <- rep(c(1:100), each=58)
boots_all_WBB$pop <- "WBB"
runs_WBB <- rbind(boots_all_WBB)
runs_WBB$plot_id <- paste(runs_WBB$pop, runs_WBB$run, sep="_")
runs_WBB$pop <- factor(runs_WBB$pop, levels = c("WBB"))
psmc_boot_WBB <- boots_all_WBB %>% group_by(run) %>% mutate(row=row_number())
psmc_boot_med_WBB <- aggregate(psmc_boot_WBB[,c("V1","V2")], list(psmc_boot_WBB$row), FUN=median)
colnames(psmc_boot_med_WBB)[1] <- "row"
psmc_boot_med_WBB$pop <- "WBB"
MED_WBB <- rbind(psmc_boot_med_WBB)
MED_WBB$pop <- factor(MED_WBB$pop, levels = c("WBB"))

########### and for NHB
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/psmc/ARRB_99_1026/boot_mu1.56")
boots_filenames_NHB <- list.files(pattern="*.txt")
boots_NHB <- lapply(boots_filenames_NHB,function(i){
  read.delim(i, header=FALSE)
})
boots_all_NHB <- do.call(rbind, boots_NHB)
boots_all_NHB$run <- rep(c(1:101), each=58)
boots_all_NHB$pop <- "NHB"
runs_NHB <- rbind(boots_all_NHB)
runs_NHB$plot_id <- paste(runs_NHB$pop, runs_NHB$run, sep="_")
runs_NHB$pop <- factor(runs_NHB$pop, levels = c("NHB"))
psmc_boot_NHB <- boots_all_NHB %>% group_by(run) %>% mutate(row=row_number())
psmc_boot_med_NHB <- aggregate(psmc_boot_NHB[,c("V1","V2")], list(psmc_boot_NHB$row), FUN=median)
colnames(psmc_boot_med_NHB)[1] <- "row"
psmc_boot_med_NHB$pop <- "WBB"
MED_NHB <- rbind(psmc_boot_med_NHB)
MED_NHB$pop <- factor(MED_NHB$pop, levels = c("NHB"))




# plot all runs & median highlighted in ggplot
# the annotate("rect"..) is adding in the LGM (last glacial maximum) and LGP (last glacial period) timeline

# set generation time (and mutation rate for axis label)
gen=21.9 

# plot
psmc_boot <-ggplot()+
  annotate("rect", xmin = 11700, xmax =115000, ymin=0, ymax=300000, alpha = 0.2, fill = "gray60")+
  geom_step(data=runs_EBB, aes(x=V1,y=V2*10000, group=plot_id), color="#1EB4C4", alpha=0.07)+
  geom_step(data=MED_EBB, aes(x=V1,y=V2*10000), color="#1EB4C4", lwd=1)+
  geom_step(data=runs_WBB, aes(x=V1,y=V2*10000, group=plot_id), color="#807DBA", alpha=0.07)+
  geom_step(data=MED_WBB, aes(x=V1,y=V2*10000), color="#807DBA", lwd=1)+
  geom_step(data=runs_NHB, aes(x=V1,y=V2*10000, group=plot_id), color="#F64F2E", alpha=0.07)+
  geom_step(data=MED_NHB, aes(x=V1,y=V2*10000), color="#F64F2E", lwd=1)+
  
  scale_x_continuous(trans='log10',breaks=c(1e3,1e4,1e5,1e6,1e7), labels=c("1e+03"=expression(10^3),"1e+04"=expression(10^4),"1e+05"=expression(10^5),"1e+06"=expression(10^6),"1e+07"=expression(10^7)))+#,limits = c(3e3,2e6))+
 # scale_y_continuous(breaks=c(0,2000,4000,6000,8000,10000,12000,14000), expand=c(0,0))+
  scale_y_continuous(trans='log10',breaks = c(1000,10000,100000),labels=c("1e+03"=expression(10^3),"1e+04"=expression(10^4), "1e+05"=expression(10^5)), limits=c(1e3,3e5))+
  annotation_logticks(sides = "lb")+
  theme_classic()+
  ylab("Effective population size")+
  xlab(paste("Years ago"))
 # geom_step(data=psmc_dat, aes(x=V1, y=V2*10000), color="red") # this line just to make sure main run lines up with bootstraps 

psmc_boot
#ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/psmc/narwhal_3groups_mu1.56e08_PSMC_boot.png", width=6,height=4,dpi=1000)



########### BOWHEAD
# Quick plot without bootstrap
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/bowhead/rightwhale_map/psmc")

psmc_dat <- read.delim("88_Pang_S51.rightwhale_2.69e-8_g52.3_plot.0.txt", header=FALSE)
psmc_dat2 <- read.delim("BM_RB_2005_001.rightwhale_2.69e-8_g52.3_plot.0.txt", header=FALSE)

# plotting just 1
ggplot() + 
  geom_step(data=psmc_dat, aes(x=V1, y=V2), col="red") + 
  geom_step(data=psmc_dat2, aes(x=V1, y=V2), col="blue") + 
  scale_x_continuous(trans='log10') + 
  annotation_logticks(sides = "lb") +
  theme_classic()


# Then load in bootstrap files (change dir first)
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/bowhead/rightwhale_map/psmc/88_Pang_boot_mu2.69")

# create list of file names to load (100 per group is too much to load individually)
boots_filenames <- list.files(pattern="*.txt")

# read data in
boots <- lapply(boots_filenames,function(i){
  read.delim(i, header=FALSE)
})

boots_all <- do.call(rbind, boots)

# every 58 rows is new boot. need column to label this as run number.
boots_all$run <- rep(c(1:101), each=58)
boots_all$pop <- "all"

# merge
runs <- rbind(boots_all)

# add unique id for label+run
runs$plot_id <- paste(runs$pop, runs$run, sep="_")

# set generation time (and mutation rate for axis label)
gen=52.3

# adjusting order to make legend by latitude
runs$pop <- factor(runs$pop, levels = c("all"))

# extract median run per group
psmc_boot <- boots_all %>% group_by(run) %>% mutate(row=row_number())
psmc_boot_med <- aggregate(psmc_boot[,c("V1","V2")], list(psmc_boot$row), FUN=median)
colnames(psmc_boot_med)[1] <- "row"
psmc_boot_med$pop <- "all"


# merge
MED <- rbind(psmc_boot_med)

# adjust label order by latitude (for legend)
MED$pop <- factor(MED$pop, levels = c("all"))


###################### second sample
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/bowhead/rightwhale_map/psmc/BM_RB_2005_001_boot_mu2.69")

# create list of file names to load (100 per group is too much to load individually)
boots_filenames2 <- list.files(pattern="*.txt")

# read data in
boots2 <- lapply(boots_filenames2,function(i){
  read.delim(i, header=FALSE)
})

boots_all2 <- do.call(rbind, boots2)

boots_all2$run <- rep(c(1:101), each=58)
boots_all2$pop <- "all"

runs2 <- rbind(boots_all2)
runs2$plot_id <- paste(runs2$pop, runs2$run, sep="_")

runs2$pop <- factor(runs2$pop, levels = c("all"))

psmc_boot2 <- boots_all2 %>% group_by(run) %>% mutate(row=row_number())
psmc_boot_med2 <- aggregate(psmc_boot2[,c("V1","V2")], list(psmc_boot2$row), FUN=median)
colnames(psmc_boot_med2)[1] <- "row"
psmc_boot_med2$pop <- "all"


# merge
MED2 <- rbind(psmc_boot_med2)

# adjust label order by latitude (for legend)
MED2$pop <- factor(MED2$pop, levels = c("all"))


# plot all runs & median highlighted in ggplot
# the annotate("rect"..) is adding in the LGM (last glacial maximum) and LGP (last glacial period) timeline

psmc_boot <-ggplot()+
   # for mu 2.69 
   annotate("rect", xmin = 11700, xmax =115000, ymin=0, ymax=70000, alpha = 0.2, fill = "gray60")+
   # for mu 1.20
  # annotate("rect", xmin = 11700, xmax =115000, ymin=0, ymax=150000, alpha = 0.2, fill = "gray60")+
  
  geom_step(data=runs2, aes(x=V1,y=V2*10000, group=plot_id), color="#a00100ff", alpha=0.07)+
  geom_step(data=MED2, aes(x=V1,y=V2*10000), color="#a00100ff", lwd=1)+
  
   geom_step(data=runs, aes(x=V1,y=V2*10000, group=plot_id), color="#be7b01", alpha=0.07)+
    geom_step(data=MED, aes(x=V1,y=V2*10000), color="#be7b01", lwd=1)+
    scale_x_continuous(trans='log10',breaks=c(1e3,1e4,1e5,1e6,1e7), labels=c("1e+03"=expression(10^3),"1e+04"=expression(10^4),"1e+05"=expression(10^5),"1e+06"=expression(10^6),"1e+07"=expression(10^7)))+#,limits = c(3e3,2e6))+
    #scale_y_continuous(breaks=c(0,2000,4000,6000,8000,10000,12000,14000), expand=c(0,0))+
  scale_y_continuous(trans='log10',breaks = c(10000,50000),labels=c("1e+04"=expression(10^4),"5e+04"=expression("5 x"~10^4)))+
    annotation_logticks(sides = "lb")+
    theme_classic()+
    ylab("Effective population size")+
    xlab(paste("Years ago"))
  # geom_step(data=psmc_dat, aes(x=V1, y=V2*10000), color="red") # this line just to make sure main run lines up with bootstraps 
  
psmc_boot

#ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/bowhead/rightwhale_map/psmc/bowhead_88_PANG_and_BM_RB_2005_001_mu2.69_PSMC_boot.png", width=6,height=4,dpi=1000)
  
