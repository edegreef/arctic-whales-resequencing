library(tidyverse)

tab <- read.table("bowhead_RWmap_allvariants.ID.removefix1.vcf.gz")

colnames(tab) <- c("CHROM", "POS","ID","REF","ALT","QUAL","FILTER", "INFO","FORMAT","88_Pang","99_01",   "ARBMGH_2002_001","BMDB_06_70","BMWG07_22","BMWG07_30","BM_01_2009","BM_CH_2000_01","BM_NSA_2008_02","BM_NSA_2009_02","BM_NSA_2009_03","BM_NSA_2010_01","BM_NSA_2011_01","BM_NSA_2011_03","BM_NSA_2012_02","BM_NSA_2012_03","BM_NSA_2014_01","BM_NSA_2020_01","BM_RB_2005_001","NSA_BM_98_01","WBF_2005_0298")


tab$`88_Pang` <- substr(tab$`88_Pang`, 0, 3)
tab$`99_01` <- substr(tab$`99_01`, 0, 3)
tab$`ARBMGH_2002_001` <- substr(tab$`ARBMGH_2002_001`, 0, 3)
tab$`BMDB_06_70` <- substr(tab$`BMDB_06_70`, 0, 3)
tab$`BMWG07_22` <- substr(tab$`BMWG07_22`, 0, 3)
tab$`BMWG07_30` <- substr(tab$`BMWG07_30`, 0, 3)
tab$`BM_01_2009` <- substr(tab$`BM_01_2009`, 0, 3)
tab$`BM_CH_2000_01` <- substr(tab$`BM_CH_2000_01`, 0, 3)
tab$`BM_NSA_2008_02` <- substr(tab$`BM_NSA_2008_02`, 0, 3)
tab$`BM_NSA_2009_02` <- substr(tab$`BM_NSA_2009_02`, 0, 3)
tab$`BM_NSA_2009_03` <- substr(tab$`BM_NSA_2009_03`, 0, 3)
tab$`BM_NSA_2010_01` <- substr(tab$`BM_NSA_2010_01`, 0, 3)
tab$`BM_NSA_2011_01` <- substr(tab$`BM_NSA_2011_01`, 0, 3)
tab$`BM_NSA_2011_03` <- substr(tab$`BM_NSA_2011_03`, 0, 3)
tab$`BM_NSA_2012_02` <- substr(tab$`BM_NSA_2012_02`, 0, 3)
tab$`BM_NSA_2012_03` <- substr(tab$`BM_NSA_2012_03`, 0, 3)
tab$`BM_NSA_2014_01` <- substr(tab$`BM_NSA_2014_01`, 0, 3)
tab$`BM_NSA_2020_01` <- substr(tab$`BM_NSA_2020_01`, 0, 3)
tab$`BM_RB_2005_001` <- substr(tab$`BM_RB_2005_001`, 0, 3)
tab$`NSA_BM_98_01` <- substr(tab$`NSA_BM_98_01`, 0, 3)
tab$`WBF_2005_0298` <- substr(tab$`WBF_2005_0298`, 0, 3)

# convert "./." to "NA"s
tab2 <- tab %>% mutate(across(where(is.character), ~na_if(., "./.")))


# create new column that checks if columns are equal, while ignoring NAs
tab3 <- tab2 %>%
  select("88_Pang","99_01", "ARBMGH_2002_001","BMDB_06_70","BMWG07_22","BMWG07_30","BM_01_2009","BM_CH_2000_01","BM_NSA_2008_02","BM_NSA_2009_02","BM_NSA_2009_03","BM_NSA_2010_01","BM_NSA_2011_01","BM_NSA_2011_03","BM_NSA_2012_02","BM_NSA_2012_03","BM_NSA_2014_01","BM_NSA_2020_01","BM_RB_2005_001","NSA_BM_98_01","WBF_2005_0298") %>%
  rowwise %>%
  #mutate(match = n_distinct(na.rm=TRUE,unlist(cur_data())) == 1) %>%
  mutate(match = n_distinct(na.rm = TRUE, c_across(everything()))) %>%
  ungroup()

# add new column to existing data frame
tab$match <- tab3$match

# make list of "TRUE" - need CHROM POS
fixed_snps <- subset(tab, match == 1)

fixed_snps_list <- fixed_snps[,1:2]
write.table(fixed_snps_list, "fixed_snps_CHROM_POS_part2_allvariants.txt", col.names =F, row.names=F, quote = FALSE)
