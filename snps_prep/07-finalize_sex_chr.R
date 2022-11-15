# Script to look at annotated bed files to create a final list of X-linked and Y-linked scaffolds determined by the M1,M2,F1,F2 combo runs. Second part of script makes upset plots to show overlap of scaffolds for the combos.

library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)

# PART 1)
# Input file names. Doing X and Y separately
bed="BOW_reference.10kb.Xlinked.bed"
type="Xlinked"
scaffold_info="BOW_reference.fasta.fai"

# Bedtools command for ref to remember order of columns
#bedtools annotate -i $ref_genome.10kb.bed -files X_M1F1.bed X_M1F2.bed X_M2F1.bed X_M2F2.bed control_M1M2.bed control_F1F2.bed > $ref_genome.10kb.Xlinked.bed

# Load data and rename scaffolds
sexlinked_windows <- read_tsv(file = bed ,col_names = F) %>% 
  rename(scaf=X1,start=X2,stop=X3,M1F1=X4,M1F2=X5, M2F1=X6, M2F2=X7, M1M2=X8, F1F2=X9)

# Want to keep scaffolds that are enriched in the 4 MF combos and not enriched in controls (not an outlier in control). Looking for non-zero values across.
# For the Y chr, do not include the FF control
filtered_windows <- sexlinked_windows %>% filter(M1F1 !=0 & M1F2 !=0 & M2F1 !=0 & M2F2 !=0 & M1M2 !=0 & F1F2 !=0)

#Y?
#filtered_windows <- sexlinked_windows %>% 
#  filter(M1F1 !=0 | M1F2 !=0 | M2F1 !=0 | M2F2 !=0) %>%
#  filter(M1M2 !=0 & F1F2 == 0)

# Add scaffold lengths
scaffold_lengths <- read_tsv(scaffold_info,col_names = c("scaf","scaffold_length"))
filtered_windows <- full_join(filtered_windows,scaffold_lengths) %>% na.omit()

# Keep only unique scaf names for final list
sex_linked_scaffolds <- filtered_windows[!duplicated(filtered_windows$scaf), ]
sex_linked_scaffolds_list <- sex_linked_scaffolds[,1]

# Check total scaffold lengths to see that it is reasonable (should be close to the checks when looking at difcover results)
chrlength <- sum(sex_linked_scaffolds$scaffold_length)
genome_length <- sum(scaffold_lengths$scaffold_length)
genome_length
chrlength/genome_length

# Save list of scaffolds to use as reference for the sex chrom in the genome!
write.table(sex_linked_scaffolds_list, paste("final_scaffolds_", type, ".txt",sep=""), col.names = FALSE, row.names = FALSE, quote=FALSE)


# PART2) - upset plots

# Load X or Y bed
sexlinked_windows <- read_tsv(file = bed ,col_names = F) %>% rename(scaf=X1,start=X2,stop=X3,M1F1=X4,M1F2=X5, M2F1=X6, M2F2=X7, M1M2=X8, F1F2=X9)

# Make dataframe into matrix with 0's and 1's for each combo
windows <- sexlinked_windows[,1:3]
sexcombo <- sexlinked_windows[,4:9]
convert01 <- sexcombo %>% mutate_if(is.numeric, ~1 * (. != 0))
table <- cbind(windows, convert01)

# Make dataframe to include by scaffold only (removing window info so no "duplicate" scaffold names). Do for each run then merge later
run1 <- table %>% filter(M1F1 !=0)
unique_run1 <- run1[!duplicated(run1$scaf), ]
unique_run1 <- select(unique_run1, scaf, M1F1)

run2 <- table %>% filter(M1F2 !=0)
unique_run2 <- run2[!duplicated(run2$scaf), ]
unique_run2 <- select(unique_run2, scaf, M1F2)

run3 <- table %>% filter(M2F1 !=0)
unique_run3 <- run3[!duplicated(run3$scaf), ]
unique_run3 <- select(unique_run3, scaf, M2F1)

run4 <- table %>% filter(M2F2 !=0)
unique_run4 <- run4[!duplicated(run4$scaf), ]
unique_run4 <- select(unique_run4, scaf, M2F2)

run5 <- table %>% filter(M1M2 !=0)
unique_run5 <- run5[!duplicated(run5$scaf), ]
unique_run5 <- select(unique_run5, scaf, M1M2)

run6 <- table %>% filter(F1F2 !=0)
unique_run6 <- run6[!duplicated(run6$scaf), ]
unique_run6 <- select(unique_run6, scaf, F1F2)

scaffold_lengths <- read_tsv(scaffold_info,col_names = c("scaf","scaffold_length"))
scaf_list <- scaffold_lengths[,1]

# Merge columns for each run
merge_scaf_list <- scaf_list %>% 
  full_join(unique_run1, by=c("scaf"="scaf")) %>% 
  full_join(unique_run2, by=c("scaf"="scaf")) %>% 
  full_join(unique_run3, by=c("scaf"="scaf")) %>% 
  full_join(unique_run4, by=c("scaf"="scaf")) %>% 
  full_join(unique_run5, by=c("scaf"="scaf")) %>% 
  full_join(unique_run6, by=c("scaf"="scaf"))

# Convert NAs to 0's
merge_scaf_list[is.na(merge_scaf_list)] <- 0
df <- as.data.frame(merge_scaf_list)

# Check that the values are correct
sum(df$M1F1)

# Convert to makecombmat format
df <- df[,-1]
df <- make_comb_mat(df, mode="intersect")

# See which combos want to extract
df[1:4]

###
###

# all 6 runs: code 111111
#df_extracted <- extract_comb(df, "111111", "111110")
#UpSet(df_extracted)

# x-linked
brewer.pal(n=6, name="Blues")

UpSet((df[comb_size(df) <= 1000]), 
      set_order=c("M1F1", "M1F2", "M2F1", "M2F2", "M1M2", "F1F2"),
      lwd=4, pt_size=unit(5, "mm"),
      comb_col=c("#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C", "black")[comb_degree((df[comb_size(df) <= 1000]))],
      border=TRUE,
      top_annotation=upset_top_annotation((df[comb_size(df) <= 1000]), annotation_name_rot=90, annotation_name_side="left", axis_param=list(side="left")),
      right_annotation=upset_right_annotation((df[comb_size(df) <= 1000]), gp=gpar(fill="gray", border=FALSE), annotation_name_side="top", axis_param=list(side="top")),
      column_title="X-linked scaffolds")


# y-linked
brewer.pal(n=5, name="Blues")
UpSet((df[comb_size(df) <= 10000]), 
      set_order=c("M1F1", "M1F2", "M2F1", "M2F2", "M1M2"),
      lwd=4, pt_size=unit(5, "mm"),
      comb_col=c("#BDD7E7", "#6BAED6", "#3182BD", "#08519C", "black")[comb_degree((df[comb_size(df) <= 10000]))],
      border=TRUE,
      top_annotation=upset_top_annotation((df[comb_size(df) <= 10000]), annotation_name_rot=90, annotation_name_side="left", axis_param=list(side="left")),
      right_annotation=upset_right_annotation((df[comb_size(df) <= 10000]), gp=gpar(fill="gray", border=FALSE), annotation_name_side="top", axis_param=list(side="top")),
      column_title="Y-linked scaffolds")
