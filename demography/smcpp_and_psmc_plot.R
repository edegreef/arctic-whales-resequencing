# Plotting SMC++ results for narwhal and bowhead, and then also plotting with PSMC results

library(tidyverse)

###################
##### SMC++ narwhal
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/smcpp")

# narwhal subgroup:
NAR_smc_EBB <- read.csv("smc_narwhal_EBB_1.56e-8_plot.csv", header=T)
NAR_smc_WBB <- read.csv("smc_narwhal_WBB_1.56e-8_plot.csv", header=T)
NAR_smc_NHB <- read.csv("smc_narwhal_RB_1.56e-8_plot.csv", header=T)

# Add unique id for label+run
NAR_smc_EBB$plot_id <- paste(NAR_smc_EBB$label, NAR_smc_EBB$plot_num, sep="_")
NAR_smc_WBB$plot_id <- paste(NAR_smc_WBB$label, NAR_smc_WBB$plot_num, sep="_")
NAR_smc_NHB$plot_id <- paste(NAR_smc_NHB$label, NAR_smc_NHB$plot_num, sep="_")

# Set generation time (and mutation rate for axis label)
NAR_gen=21.9

NAR_EBB <- NAR_smc_EBB %>% group_by(plot_num) %>% mutate(row=row_number())
NAR_WBB <- NAR_smc_WBB %>% group_by(plot_num) %>% mutate(row=row_number())
NAR_NHB <- NAR_smc_NHB %>% group_by(plot_num) %>% mutate(row=row_number())

NAR_med_EBB <- aggregate(NAR_EBB[,c("x","y")], list(NAR_EBB$row), FUN=median)
colnames(NAR_med_EBB)[1] <- "row"
NAR_med_EBB$label <- "EBB"

NAR_med_WBB <- aggregate(NAR_WBB[,c("x","y")], list(NAR_WBB$row), FUN=median)
colnames(NAR_med_WBB)[1] <- "row"
NAR_med_WBB$label <- "WBB"

NAR_med_NHB <- aggregate(NAR_NHB[,c("x","y")], list(NAR_NHB$row), FUN=median)
colnames(NAR_med_NHB)[1] <- "row"
NAR_med_NHB$label <- "NHB"


# plot iterations and highlight median
ggplot()+
  annotate("rect", xmin = 11700, xmax =115000, ymin=2000, ymax=50000, alpha = 0.5, fill = "gray88")+
  #geom_line(data=NAR_smc, aes(x=(x*NAR_gen),y=y, group=plot_id, alpha=0.5), colour="#c8e1f7")+
  geom_line(data=NAR_med_EBB, aes(x=(x*NAR_gen),y=y), lwd=1.5, alpha=1, colour="#6baed6")+
  geom_line(data=NAR_med_WBB, aes(x=(x*NAR_gen),y=y), lwd=1.5, alpha=1, colour="#2171B5")+
  geom_line(data=NAR_med_NHB, aes(x=(x*NAR_gen),y=y), lwd=1.5, alpha=1, colour="#bdd7e7")+
  
  xlab("Years ago")+
  ylab("Effective population size")+
  annotation_logticks(sides = "lb")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(trans='log10', breaks=c(10000,100000), limits=c(2000,50000), labels=c("1e+04"=expression(10^4), "1e+05"=expression(10^5)), expand=c(0,0))+
  # ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  scale_x_continuous(breaks=c(1e1,1e2,1e3,1e4,1e5,1e6),trans='log10',labels=c("1e+01"=expression(10^1),"1e+02"=expression(10^2),"1e+03"=expression(10^3), "1e+04"=expression(10^4), "1e+05"=expression(10^5), "1e+06"=expression(10^6)), limits=c(1000,4000000),expand=c(0,0))+
  theme(legend.position = "none")

###################
##### PSMC narwhal
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/narwhal/psmc/")

NAR_psmc_EBB <- read.delim("ARCR_07_1065/ARCR_07_1065_1.56e-8_g21.9plot.0.txt", header=FALSE)
NAR_psmc_EBB$ne <- NAR_psmc_EBB$V2 * 10000
# b/c psmc not great at recent stuff, maybe remove the most recent 5 rows?? just want to show the deeper past anyway 
NAR_psmc_EBB2 <- NAR_psmc_EBB[7:58,]

NAR_psmc_WBB <- read.delim("ARGF_07_1127/ARGF_07_1127_1.56e-8_g21.9plot.0.txt", header=FALSE)
NAR_psmc_WBB$ne <- NAR_psmc_WBB$V2 * 10000
NAR_psmc_WBB2 <- NAR_psmc_WBB[7:58,]

NAR_psmc_NHB <- read.delim("ARRB_99_1026/ARRB_99_1026_1.56e-8_g21.9plot.0.txt", header=FALSE)
NAR_psmc_NHB$ne <- NAR_psmc_NHB$V2 * 10000
NAR_psmc_NHB2 <- NAR_psmc_NHB[7:58,]

ggplot()+
  # ice age
    annotate("rect", xmin = 11700, xmax =2500000, ymin=2000, ymax=50000, alpha = 0.5, fill = "gray92")+
  # LGP
  annotate("rect", xmin = 11700, xmax =115000, ymin=2000, ymax=50000, alpha = 0.5, fill = "gray85")+

  #stat_smooth(data=NAR_med, aes(x=(x*NAR_gen),y=y), method=lm, formula = y ~ poly(x,10), se=FALSE,  colour="#2171B5", lwd=1.2)+ 
  #geom_line(data=NAR_med, aes(x=(x*NAR_gen),y=y), lwd=1.5, alpha=1, colour="orange")+
  
  stat_smooth(data=NAR_psmc_EBB2, aes(x=V1,y=ne),method=lm, formula=y~poly(x,10),se=FALSE, lwd=1.2, colour="#6baed6", linetype="dashed")+
  stat_smooth(data=NAR_psmc_WBB2, aes(x=V1,y=ne),method=lm, formula=y~poly(x,10),se=FALSE, lwd=1.2, colour="#2171B5", linetype="dashed")+
  stat_smooth(data=NAR_psmc_NHB2, aes(x=V1,y=ne),method=lm, formula=y~poly(x,10),se=FALSE, lwd=1.2, colour="#bdd7e7", linetype="dashed")+
  
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


###################
##### SMC++ bowhead

setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/")

BOW_smc <- read.csv("bowhead/rightwhale_map/smcpp/smc_RWmap_2.69e-8_plot.csv", header=T)

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
  geom_line(data=BOW_smc, aes(x=(x*BOW_gen),y=y, group=plot_id, alpha=0.5), colour="orange")+
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


###################
##### PSMC bowhead

BOW_psmc <- read.delim("bowhead/rightwhale_map/psmc/88_Pang_S51.rightwhale_2.69e-8_g52.3_plot.0.txt", header=FALSE)
BOW_psmc$ne <- BOW_psmc$V2 * 10000

BOW_psmc2 <- BOW_psmc[7:58,]

###################
######### Combine SMC++ and PSMC for narwhal and bowhead
ggplot()+
  # ice age
  annotate("rect", xmin = 11700, xmax =2500000, ymin=2000, ymax=60000, alpha = 0.5, fill = "gray92")+
  # LGP
  annotate("rect", xmin = 11700, xmax =115000, ymin=2000, ymax=60000, alpha = 0.5, fill = "gray85")+
  stat_smooth(data=BOW_med, aes(x=(x*BOW_gen),y=y), method=lm, formula = y ~ poly(x,7), se=FALSE,  colour="#F8D568", lwd=1.2)+ 
  #geom_line(data=BOW_psmc2, aes(x=V1,y=ne), lwd=1.5, alpha=1, colour="yellow")+
  #geom_line(data=BOW_med, aes(x=(x*BOW_gen),y=y), lwd=1.5, alpha=1, colour="yellow")+
  
  
  stat_smooth(data=BOW_psmc2, aes(x=V1,y=ne),method=lm, formula=y~poly(x,10),se=FALSE, lwd=1.2, colour="#D59521", linetype="dashed")+
  
  ### added in narwhal here too

   stat_smooth(data=NAR_med_WBB, aes(x=(x*NAR_gen),y=y), method=lm, formula = y ~ poly(x,10), se=FALSE,  colour="#2171B5", lwd=1.2)+ 
   stat_smooth(data=NAR_med_NHB, aes(x=(x*NAR_gen),y=y), method=lm, formula = y ~ poly(x,10), se=FALSE,  colour="#bdd7e7", lwd=1.2)+ 
   stat_smooth(data=NAR_med_EBB, aes(x=(x*NAR_gen),y=y), method=lm, formula = y ~ poly(x,10), se=FALSE,  colour="#6baed6", lwd=1.2)+ 
   stat_smooth(data=NAR_psmc_WBB2, aes(x=V1,y=ne),method=lm, formula=y~poly(x,10),se=FALSE, lwd=1.2, colour="#2171B5", linetype="dashed")+
   stat_smooth(data=NAR_psmc_NHB2, aes(x=V1,y=ne),method=lm, formula=y~poly(x,10),se=FALSE, lwd=1.2, colour="#bdd7e7", linetype="dashed")+
   stat_smooth(data=NAR_psmc_EBB2, aes(x=V1,y=ne),method=lm, formula=y~poly(x,10),se=FALSE, lwd=1.2, colour="#6baed6", linetype="dashed")+
  ### back to plot
  xlab("Thousand years ago (kya)")+
  ylab(expression(paste("Effective population size (")*10^{3}*")"))+
  annotation_logticks(sides = "lb")+
  theme_bw()+
  #scale_y_continuous(trans='log10', breaks=c(10000,100000), limits=c(2000,100000), labels=c("1e+04"=expression(10^4), "1e+05"=expression(10^5)), expand=c(0,0))+
  # ylab(expression(paste("Effective population size (")~10^{3}~")"))+
  scale_x_continuous(breaks=c(1e1,1e2,1e3,1e4,1e5,1e6,1e7),trans='log10',labels=c("1e+01"=expression(10^1),"1e+02"=expression(10^2),"1e+03"=expression(1), "1e+04"=expression(10), "1e+05"=expression(100), "1e+06"="1,000","1e+07"="10,000"), limits=c(1000,20000000),expand=c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous(trans='log10', breaks=c(10000,100000), limits=c(2000,60000), labels=c("1e+04"=expression(10), "1e+05"=expression(100)), expand=c(0,0))+
  theme(legend.position = "none")

#ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/demo_history/plots/SMCPP_PSMC_narwhal_subgroups_bowhead_RWmap.png", height=4, width=5.5, dpi=1000)

