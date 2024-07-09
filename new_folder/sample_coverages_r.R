setwd("C:/Users/eveli/Dropbox/Whales with Garroway/02-killerwhale")

library(openxlsx)
library(ggplot2)

# load sample info 
sample_info <- read.xlsx("Garroway-DeGreef UofM killerwhale_genomic_sample_info_23Sept19.xlsx")

sample_info <- subset(sample_info, genome_sample_ID != "MM406")
sample_info$row <- 1:nrow(sample_info)

# modal coverage
bar_mod <- ggplot(data=sample_info, aes(x=factor(row, labels=genome_sample_ID), y=modal_coverage))+
  geom_bar(stat="identity", fill="gray")+
  theme_classic()+
  xlab("individual")+
  ylab("modal coverage")+
 # ggtitle("individual modal coverage")+
  theme(plot.title=element_text(hjust=0.5),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1))+
  geom_hline(yintercept=19, color="red")

  #coord_cartesian(expand=FALSE)
#scale_fill_manual( values = c( "no"="gray", "yes"="orange" ), guide = FALSE )
bar_mod


# mean coverage
bar_mean <- ggplot(data=sample_info, aes(x=factor(row, labels=genome_sample_ID), y=mean_coverage))+
  geom_bar(stat="identity", fill="gray")+
  theme_classic()+
  xlab("individual")+
  ylab("mean coverage")+
  # ggtitle("individual modal coverage")+
  theme(plot.title=element_text(hjust=0.5),
        #axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1))
#coord_cartesian(expand=FALSE)
#scale_fill_manual( values = c( "no"="gray", "yes"="orange" ), guide = FALSE )
bar_mean

library(patchwork)
bar_mod / bar_mean

mean(sample_info$modal_coverage)
mean(sample_info$mean_coverage)
