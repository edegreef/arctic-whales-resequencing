# Plot epos .dat files

library(tidyverse)
library(ggplot2)
library(patchwork)

# Here using Narwhal example

# Load .dat file
NAR_all <- read_delim("NAR/snps/epos/NAR_frq_epos_boot10000_mu.dat", col_names=c("Gen",	"LowerQ",	"Median",	"UpperQ"), skip = 1)

# Also going to plot an epos run using a subset of samples (22/58) to see how it varies
NAR_22 <- read_delim("NAR/snps/epos/NAR_frq_epos_boot_mu_n22.dat", col_names=c("Gen",	"LowerQ",	"Median",	"UpperQ"), skip = 1)

# Quick plot
ggplot()+
  geom_line(data=NAR, aes(x=Gen, y=LowerQ), colour="purple", lwd=1)+
  geom_line(data=NAR, aes(x=Gen, y=UpperQ), colour="lightblue", lwd=1)+
  geom_line(data=NAR, aes(x=Gen, y=Median), lwd=1.5)+
  theme_bw()
  
# Remember "time" in the epos file is "generations"
# Generation times:
# 21.9 for NAR 
# 52.3 for BOW

NAR_all$time <- (NAR_all$Gen * 21.9)
NAR_22$time <- (NAR_22$Gen * 21.9)

NAR_all_plot <- ggplot()+
  geom_ribbon(data=NAR_all, aes(x=time, ymin = LowerQ, ymax = UpperQ), colour = "gray", alpha = 0.2) +
  geom_line(data=NAR_all, aes(x=time, y=Median), lwd=1)+ 
  geom_line(data=NAR_all, aes(x=time, y=LowerQ), colour="purple", lwd=.7, linetype="dashed")+
  geom_line(data=NAR_all, aes(x=time, y=UpperQ), colour="lightblue", lwd=.7, linetype="dashed")+
  theme_bw()+
  scale_x_continuous(breaks=c(1e4,1e5,1e6,1e7,1e8,1e9),trans='log10',limits=c(1e3,3e8))+
  scale_y_continuous(trans='log10', limits=c(1e4,5e7), expand=c(0,0))+
  annotation_logticks(sides = "lb")+
  xlab("Years ago")+
  ylab("Effective Pop Size")+
  ggtitle("Narwhal n=58")

#ggsave("BOW_epos_10000boot_log.png",width=5,height=4,dpi=400)

NAR_22_plot <- ggplot()+
  geom_ribbon(data=NAR_22, aes(x=time, ymin = LowerQ, ymax = UpperQ), colour = "gray", alpha = 0.2) +
  geom_line(data=NAR_22, aes(x=time, y=Median), lwd=1)+ 
  geom_line(data=NAR_22, aes(x=time, y=LowerQ), colour="purple", lwd=.7, linetype="dashed")+
  geom_line(data=NAR_22, aes(x=time, y=UpperQ), colour="lightblue", lwd=.7, linetype="dashed")+
  theme_bw()+
  scale_x_continuous(breaks=c(1e4,1e5,1e6,1e7,1e8,1e9),trans='log10',limits=c(1e3,3e8))+
  scale_y_continuous(trans='log10', limits=c(1e4,5e7), expand=c(0,0))+
  annotation_logticks(sides = "lb")+
  xlab("Years ago")+
  ylab("Effective Pop Size")+
  ggtitle("Narwhal subset n=22")

NAR_all_plot / NAR_22_plot
#ggsave("NAR_epos_1000boot_log_ABonly.png",width=5,height=4,dpi=400)


