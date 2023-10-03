# Plot PSMC data with and without bootstrapping

library(tidyverse)

# First do without bootstrap
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/bowhead/psmc")

psmc_singlerun <- read.delim("88_Pang_2.69e-8_g52.3_plot.0.txt", header=FALSE)

# plotting just 1
ggplot(data=psmc_singlerun, aes(x=V1, y=V2)) + 
  geom_step() + 
  scale_x_continuous(breaks=c(1e3,1e4,1e5,1e6),trans='log10',limits=c(1e4,1e6)) + 
  scale_y_continuous(limits=c(0,2))+ 
  #geom_hline(yintercept = 5) + 
  annotation_logticks(sides = "lb") +
  theme_classic()


# With boostrapping files
# create list of file names to load (100 per group is too much to load individually)
psmc_filenames <- list.files(pattern="combined*")

psmc_files <- lapply(psmc_filenames,function(i){
  read.delim(i, header=FALSE)
})

psmc_boot <- do.call(rbind, psmc_files)

# Every 58 rows is new boot. need column to label this as run number.
psmc_boot$run <- rep(c(1:101), each=58)

psmc_boot$pop <- "all"

# Merge
runs <- rbind(psmc_boot)

# Add unique id for label+run
runs$plot_id <- paste(runs$pop, runs$run, sep="_")

# Set generation time (21.9 for narwhal, 52.3 for bowhead)
gen=52.3

# adjusting order to make legend by latitude
runs$pop <- factor(runs$pop, levels = c("all"))

# extract median run per group
psmc <- psmc_boot %>% group_by(run) %>% mutate(row=row_number())
psmc_med <- aggregate(psmc[,c("V1","V2")], list(psmc$row), FUN=median)
colnames(psmc_med)[1] <- "row"
psmc_med$pop <- "all"

# Merge
MED <- rbind(psmc_med)

# Adjust label order by latitude (for legend)
MED$pop <- factor(MED$pop, levels = c("all"))

# Plot all runs & median highlighted in ggplot
psmc_boot <-ggplot()+
  # Highlight LGP
  annotate("rect", xmin = 11700, xmax =115000, ymin=0, ymax=40000, alpha = 0.2, fill = "gray60")+
  # Plot bootstraps
  geom_step(data=runs, aes(x=V1,y=V2*10000, group=plot_id), color="#be7b01", alpha=0.1)+
  # Plot median run
  geom_step(data=MED, aes(x=V1,y=V2*10000), color="#be7b01", lwd=1.5)+
  #scale_color_manual(values=group_colors)+
  scale_x_continuous(trans='log10',breaks=c(1e3,1e4,1e5,1e6,1e7), labels=c("1e+03"=expression(10^3),"1e+04"=expression(10^4),"1e+05"=expression(10^5),"1e+06"=expression(10^6),"1e+07"=expression(10^7)))+#,limits = c(3e3,2e6))+
  #scale_y_continuous(breaks=c(0,2000,4000,6000,8000,10000,12000,14000), expand=c(0,0))+
  scale_y_continuous(trans='log10')+# limits = c(0.1,50))+
  annotation_logticks(sides = "lb")+
  theme_classic()+
  ylab("Effective population size")+
  xlab(paste("Years ago"))

psmc_boot
#ggsave("bowhead_PSMC_boot.png", width=6,height=4,dpi=1000)
