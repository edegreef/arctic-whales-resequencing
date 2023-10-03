# Plotting snmf/admixture results into ggplot

# Load modules
library(tidyverse)
library(reshape2)
library(patchwork)

# Directory
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/narwhal/LEA/pop structure/location_corrected/")

# Load qmatrix (from snmf_lea step earlier)
qmatrix <- read.csv("qmatrix_K2_poporder.csv")
qmatrix <- subset(qmatrix, select = -c(pop_order))

# Convert dataframe to long format
qlong = melt(qmatrix, id.vars=c("Ind","Site"))
head(qlong)

snmf_admix = ggplot(data=qlong, aes(x=fct_inorder(Ind), y=value, fill=variable))+
  geom_bar(stat = "identity",position = position_stack(reverse = TRUE), width=1)+
  scale_y_continuous(expand = c(0,0))+
  ylab("K=2")+
  xlab("")+
  theme_classic()+
  scale_fill_manual(values = c("Cluster.1" = "#F64F2E", "Cluster.2" = "#807DBA"))+
 # theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank())+theme(legend.position = "none")+
  coord_cartesian(expand=FALSE)+
  theme(axis.text.x=element_blank())
# theme(plot.title=element_text(hjust=0.5, face="bold"),axis.text.x=element_text(angle=90))

snmf_admix

## Now adding K=3 and K=5
qmatrix_K3 <- read.csv("qmatrix_K3_poporder.csv")
qmatrix_K3 <- subset(qmatrix_K3, select = -c(pop_order))

#qmatrix_K4 <- read.csv("qmatrix_K4_poporder.csv")
#qmatrix_K4 <- subset(qmatrix_K4, select = -c(pop_order))

# convert dataframe to long format
qlong_K3 = melt(qmatrix_K3, id.vars=c("Ind","Site"))
head(qlong_K3)


## Plots
snmf_admix_K3 = ggplot(data=qlong_K3, aes(x=fct_inorder(Ind), y=value, fill=variable))+
  geom_bar(stat = "identity",position = position_stack(reverse = TRUE), width=1)+
  scale_y_continuous(expand = c(0,0))+
  ylab("K=3")+
  xlab("")+
  theme_classic()+
  scale_fill_manual(values = c("Cluster.1" = "#F64F2E", "Cluster.3" = "#807DBA", "Cluster.2" = "#1EB4C4"))+
  # theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank())+theme(legend.position = "none")+
  coord_cartesian(expand=FALSE)+
  theme(axis.text.x=element_blank())+
  theme(plot.title=element_text(hjust=0.5, face="bold"),axis.text.x=element_text(angle=90))
snmf_admix_K3

#snmf_admix_K4 = ggplot(data=qlong_K4, aes(x=fct_inorder(Ind), y=value, fill=variable))+
#  geom_bar(stat = "identity",position = position_stack(reverse = TRUE), width=1)+
#  scale_y_continuous(expand = c(0,0))+
#  ylab("K=4")+
#  xlab("Individual")+
#  theme_classic()+
#  scale_fill_manual(values = c("Cluster.3" = "#F64F2E", "Cluster.2" = "#807DBA", "Cluster.1" = #"#1EB4C4", "Cluster.4" = "#FED976"))+
  # theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank())+theme(legend.position = "none")+
 # coord_cartesian(expand=FALSE)+
  #theme(axis.text.x=element_blank())+
#  theme(plot.title=element_text(hjust=0.5, face="bold"),axis.text.x=element_text(angle=90))
#snmf_admix_K4

snmf_admix/snmf_admix_K3

ggsave("K2-3_barplots_loccorrected.svg", width=7, height=4, dpi=1000)
