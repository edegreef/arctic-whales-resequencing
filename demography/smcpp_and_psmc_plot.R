# Plot bowhead and narwhal separately for DFO report
# Plot bowhead and narwhal together for paper

library(tidyverse)


setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/smcpp")

###########
# NARWHAL:
NAR_smc <- read.csv("smc_narwhal_1group_100iter_2DL_plot.csv", header=T)

# Add unique id for label+run
NAR_smc$plot_id <- paste(NAR_smc$label, NAR_smc$plot_num, sep="_")

# Set generation time (and mutation rate for axis label)
NAR_gen=21.9
#NAR_mu=paste("1.56x",expression(10^{-8}), sep="")

NAR <- NAR_smc %>% 
  group_by(plot_num) %>% 
  mutate(row=row_number())
NAR_med <- aggregate(NAR[,c("x","y")], list(NAR$row), FUN=median)
colnames(NAR_med)[1] <- "row"
NAR_med$label <- "NAR"


# Plot iterations and highlight median
ggplot()+
  annotate("rect", xmin = 11700, xmax =115000, ymin=2000, ymax=50000, alpha = 0.5, fill = "gray88")+
  geom_line(data=NAR_smc, aes(x=(x*NAR_gen),y=y, group=plot_id, alpha=0.5), colour="#c8e1f7")+
  geom_line(data=NAR_med, aes(x=(x*NAR_gen),y=y), lwd=1.5, alpha=1, colour="#2171B5")+
  xlab("Years ago")+
  ylab("Effective population size")+
  annotation_logticks(sides = "lb")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(trans='log10', breaks=c(10000,100000), limits=c(2000,50000), labels=c("1e+04"=expression(10^4), "1e+05"=expression(10^5)), expand=c(0,0))+
  # ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  scale_x_continuous(breaks=c(1e1,1e2,1e3,1e4,1e5,1e6),trans='log10',labels=c("1e+01"=expression(10^1),"1e+02"=expression(10^2),"1e+03"=expression(10^3), "1e+04"=expression(10^4), "1e+05"=expression(10^5), "1e+06"=expression(10^6)), limits=c(1000,4000000),expand=c(0,0))+
  theme(legend.position = "none")

###########
# Bring in PSMC
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/psmc")
# psmc
NAR_psmc <- read.delim("ARCR_07_1065_1.56e-8_g21.9_plot.0.txt", header=FALSE)
NAR_psmc$ne <- NAR_psmc$V2 * 10000

# quick plot for psmc (not log scaled) ### remember PSMC already has gen and mu taken into account
ggplot()+
  geom_line(data=NAR_psmc, aes(x=V1,y=ne), lwd=1.2, colour="royalblue")


# b/c psmc not great at recent stuff, maybe remove the most recent 5 rows?? just want to show the deeper past anyway 
NAR_psmc2 <- NAR_psmc[6:58,]

ggplot()+
  
  # ice age
    annotate("rect", xmin = 11700, xmax =2500000, ymin=2000, ymax=50000, alpha = 0.5, fill = "gray92")+
  # LGP
  annotate("rect", xmin = 11700, xmax =115000, ymin=2000, ymax=50000, alpha = 0.5, fill = "gray85")+

  #stat_smooth(data=NAR_med, aes(x=(x*NAR_gen),y=y), method=lm, formula = y ~ poly(x,10), se=FALSE,  colour="#2171B5", lwd=1.2)+ 
  #geom_line(data=NAR_med, aes(x=(x*NAR_gen),y=y), lwd=1.5, alpha=1, colour="orange")+
  stat_smooth(data=NAR_psmc2, aes(x=V1,y=ne),method=lm, formula=y~poly(x,10),se=FALSE, lwd=1.2, colour="#08306B", linetype="dashed")+
  xlab("Thousand years ago (kya)")+
  ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  annotation_logticks(sides = "lb")+
  theme_bw()+
  #scale_y_continuous(trans='log10', breaks=c(10000,100000), limits=c(2000,100000), labels=c("1e+04"=expression(10^4), "1e+05"=expression(10^5)), expand=c(0,0))+
  # ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  scale_x_continuous(breaks=c(1e1,1e2,1e3,1e4,1e5,1e6,1e7),trans='log10',labels=c("1e+01"=expression(10^1),"1e+02"=expression(10^2),"1e+03"=expression(1), "1e+04"=expression(10), "1e+05"=expression(100), "1e+06"="1,000","1e+07"="10,000"), limits=c(1000,5000000),expand=c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(trans='log10', breaks=c(10000,100000), limits=c(2000,50000), labels=c("1e+04"=expression(10), "1e+05"=expression(100)), expand=c(0,0))+
  theme(legend.position = "none")
# ggtitle("SMC++")

#ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/demo_history/plots/SMCPP_PSMC_narwhal_logged_v2.2.png", height=3, width=4, dpi=1000)
#####


######## 
# BOWHEAD WHALE:
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/bowhead/smcpp")

BOW_smc <- read.csv("smc_bowhead_1group_100iter_2DL_plot.csv", header=T)

# Add unique id for label+run
BOW_smc$plot_id <- paste(BOW_smc$label, BOW_smc$plot_num, sep="_")

# Set generation time (and mutation rate for axis label)
BOW_gen=52.3
#NAR_mu=paste("1.56x",expression(10^{-8}), sep="")

BOW <- BOW_smc %>% 
  group_by(plot_num) %>% 
  mutate(row=row_number())
BOW_med <- aggregate(BOW[,c("x","y")], list(BOW$row), FUN=median)
colnames(BOW_med)[1] <- "row"
BOW_med$label <- "BOW"


# plot iterations and highlight median
ggplot()+
  annotate("rect", xmin = 11700, xmax =115000, ymin=2000, ymax=50000, alpha = 0.5, fill = "gray88")+
  #geom_line(data=NAR_smc, aes(x=(x*NAR_gen),y=y, group=plot_id, alpha=0.5), colour="#c8e1f7")+
  geom_line(data=BOW_med, aes(x=(x*BOW_gen),y=y), lwd=1.5, alpha=1, colour="orange")+
  xlab("Years ago")+
  ylab("Effective population size")+
  annotation_logticks(sides = "lb")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(trans='log10', breaks=c(10000,100000), limits=c(2000,50000), labels=c("1e+04"=expression(10^4), "1e+05"=expression(10^5)), expand=c(0,0))+
  # ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  scale_x_continuous(breaks=c(1e1,1e2,1e3,1e4,1e5,1e6),trans='log10',labels=c("1e+01"=expression(10^1),"1e+02"=expression(10^2),"1e+03"=expression(10^3), "1e+04"=expression(10^4), "1e+05"=expression(10^5), "1e+06"=expression(10^6)), limits=c(1000,4000000),expand=c(0,0))+
  theme(legend.position = "none")

###########
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/bowhead/psmc")
# psmc
# psmc
BOW_psmc <- read.delim("88_Pang_2.69e-8_g52.3_plot.0.txt", header=FALSE)
BOW_psmc$ne <- BOW_psmc$V2 * 10000

BOW_psmc2 <- BOW_psmc[6:58,]

# quick plot for psmc (not log scaled)
ggplot()+
  geom_line(data=BOW_psmc2, aes(x=V1,y=ne), lwd=1.2, colour="#be7b01")


## Plotting SMC++ and PSMC for both narwhal and bowhead whale together:
ggplot()+
  # Ice age
  annotate("rect", xmin = 11700, xmax =2500000, ymin=2000, ymax=50000, alpha = 0.5, fill = "gray92")+
  # LGP
  annotate("rect", xmin = 11700, xmax =115000, ymin=2000, ymax=50000, alpha = 0.5, fill = "gray85")+
  stat_smooth(data=BOW_med, aes(x=(x*BOW_gen),y=y), method=lm, formula = y ~ poly(x,7), se=FALSE,  colour="orange", lwd=1.2)+ 
 # geom_line(data=BOW_med, aes(x=(x*BOW_gen),y=y), lwd=1.5, alpha=1, colour="yellow")+
  stat_smooth(data=BOW_psmc2, aes(x=V1,y=ne),method=lm, formula=y~poly(x,10),se=FALSE, lwd=1.2, colour="#be7b01", linetype="dashed")+
  ### added in narwhal here too
  stat_smooth(data=NAR_med, aes(x=(x*NAR_gen),y=y), method=lm, formula = y ~ poly(x,10), se=FALSE,  colour="#2171B5", lwd=1.2)+ 
  #geom_line(data=NAR_med, aes(x=(x*NAR_gen),y=y), lwd=1.5, alpha=1, colour="orange")+
  stat_smooth(data=NAR_psmc2, aes(x=V1,y=ne),method=lm, formula=y~poly(x,10),se=FALSE, lwd=1.2, colour="#08306B", linetype="dashed")+
  
  xlab("Thousand years ago (kya)")+
  ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  annotation_logticks(sides = "lb")+
  theme_bw()+
  #scale_y_continuous(trans='log10', breaks=c(10000,100000), limits=c(2000,100000), labels=c("1e+04"=expression(10^4), "1e+05"=expression(10^5)), expand=c(0,0))+
  # ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  scale_x_continuous(breaks=c(1e1,1e2,1e3,1e4,1e5,1e6,1e7),trans='log10',labels=c("1e+01"=expression(10^1),"1e+02"=expression(10^2),"1e+03"=expression(1), "1e+04"=expression(10), "1e+05"=expression(100), "1e+06"="1,000","1e+07"="10,000"), limits=c(1000,20000000),expand=c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(trans='log10', breaks=c(10000,100000), limits=c(2000,50000), labels=c("1e+04"=expression(10), "1e+05"=expression(100)), expand=c(0,0))+
  theme(legend.position = "none")
 # theme(text=element_text(family="serif"))


# original size:
#ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/plots/SMCPP_PSMC_bowhaed_and_narwhal_logged.png", height=3, width=4, dpi=1000)

ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/plots/SMCPP_PSMC_bowhaed_and_narwhal_logged_size2_serif.png", height=4, width=5.5, dpi=1000)


