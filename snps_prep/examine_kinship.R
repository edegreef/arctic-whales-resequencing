# looking at kinship estimates

library(plinkQC)
library(ggplot2)

directory="C:/Users/Evelien de Greef/Dropbox/Whales with Garroway/BOW/snps/kinship"
files="BOW_SNPS.filter1.miss.biallel.min100kb.autosomes.maf05.LDprunedr08" ##for .genome and .imiss
  
  
# Quick plot
#png("plinkQC_BOW.png", w=1000, h=500, res=100)
evaluate_check_relatedness(
  qcdir=directory,
  name=files,
  highIBDTh = 0.1875,
  imissTh = 0.03,
  interactive = FALSE,
  verbose = FALSE
)
#dev.off()
 
# Extract values for the related individuals from the .genome file (filtering by pi_hat values)
IBD <- read.table(paste(directory,"/", files,".genome", sep=""), header=T)
IBD_order_PH <- IBD[order(IBD[,10], decreasing=TRUE),]
write.csv(IBD_order_PH, "IBD_BOW_PIHAT.csv")

# Make nicer histogram

# Order by pi_hat and then add nrow (for pair#s)
IBD_2 <- IBD[order(IBD[,10], decreasing=TRUE),]
IBD_2$pair <- 1:nrow(IBD_2)

ggplot(data=IBD_2, aes(x=PI_HAT))+
  geom_histogram(bins=125, color="steelblue", fill="steelblue")+
  theme_bw()+
  ylab("number of pairs")+
  xlab("pi hat")
ggsave("IBD_BOW_geomhistogram.png", width=6,height=4,dpi=300)


# PI_HAT value reference:
# 1st degree relative: 0.5
# 2nd degree relative: 0.25
# 3rd degree relative: 0.125

