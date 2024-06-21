# Calculating proportion of observed heterozygosity using the .het outputs from vcftools
# vcftools --vcf input.vcf --het --out output.het

library(tidyverse)

setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/het")

# Bowhead
het_bow <- read.table("bowhead_RWmap.n20.output_het.het", header=T)

het_bow$e.het <- (het_bow$N_SITES- het_bow$E.HOM.)
het_bow$o.het <- (het_bow$N_SITES- het_bow$O.HOM.)

het_bow$e.het.prop <- (het_bow$e.het/het_bow$N_SITES)
het_bow$o.het.prop <- (het_bow$o.het/het_bow$N_SITES)

het_bow$species <- "Bowhead whale"

# Narwhal
het_nar <- read.table("narwhal.n57.output_het.het", header=T)

het_nar$e.het <- (het_nar$N_SITES- het_nar$E.HOM.)
het_nar$o.het <- (het_nar$N_SITES- het_nar$O.HOM.)

het_nar$e.het.prop <- (het_nar$e.het/het_nar$N_SITES)
het_nar$o.het.prop <- (het_nar$o.het/het_nar$N_SITES)

het_nar$species <- "Narwhal"

########

# combine
het_all <- rbind(het_bow, het_nar)

# plot
ggplot(het_all, aes(x=(species), y = o.het.prop)) + 
  geom_boxplot(width=0.5, lwd=1, fill=NA, colour="gray70")+
  geom_jitter(width=0.15, aes(colour=species), alpha=0.8,)+
  theme_bw()+
  xlab("")+
  ylab("Observed heterozygosity")+
  scale_colour_manual(name="Species",values = c("Narwhal"="#4292C6", "Bowhead whale"="#f8d568"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(), legend.position="none", axis.text.x=element_text(angle=0))+
  scale_x_discrete(limits=c("Narwhal", "Bowhead whale"))

ggsave("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/het/het_plot_updated_RWmap.png", height=3, width=3, dpi=400)
