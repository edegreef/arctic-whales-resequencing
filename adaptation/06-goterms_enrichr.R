# Looking at Magma results- prepping gene list and then using enrichR for go-term analysis
# Help from Matt Thorstensen for initial set up

library(enrichR)
library(tidyverse)

setwd("C:/Users/eveli/Dropbox/Whales with Garroway/01-arctic_whales/magma")

## Bowhead whale
# Read in magma outputs 
bow_annot <- read_delim("bowhead_magma_annotation.genes.annot", delim = "\t", skip = 3, col_names = c("geneID", "gene_loc", "snp", "snp2", "snp3", "snp4"))

# For bowhead, the gene names will already be in gene ID
# Re-format gene list
split <- as.data.frame(str_split_fixed(bow_annot$geneID, "-", 3))
colnames(split) <- c("scaffold", "hitnumber", "gene")

# Unique hits
bow_unique_info <- split[!duplicated(split$gene), ]
# Save this one for table
write.csv(bow_unique_info, "bow_unique_info.csv")

# Just gene ID
bow_unique <- as.data.frame(unique(split$gene))
colnames(bow_unique) <- "gene"


## Narwhal
# Read in magma outputs 
nar_annot <- read_delim("narwhal_magma_annotation.genes.annot", delim = "\t", skip = 3, col_names = c("geneID", "gene_loc", "snp", "snp2", "snp3", "snp4"))

# Read in gff
gff <- read_delim("C:/Users/eveli/Dropbox/Whales with Garroway/ref_genomes/annotations/Monodon_monoceros.NGI_Narwhal_1.107.gff3", delim="\t",comment="#", col_names = F)
colnames(gff) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
gff_genes_only <- subset(gff, feature=="gene")

# Want to make a column of just the geneID to match the annot file, then can left_join or something
# Pulling out the ENSMMNG
gff_genes_only[c("gene_ID", "gene_info", "gene_ID2", "logic_name", "version")] <- str_split_fixed(gff_genes_only$attribute, ';',5)
gff_genes_only[c("ID", "ENSMMNG")] <- str_split_fixed(gff_genes_only$gene_ID, ':',2)

# Since only 31 gene IDs are needed, let's try this way
nar_annot_genes <- as.data.frame(nar_annot$geneID)
colnames(nar_annot_genes) <- "geneID"

# left_join time
nar_annot_results <- left_join(nar_annot_genes, gff_genes_only, by=c("geneID"="ENSMMNG"))

# subset
nar_annot_results <- nar_annot_results[c("geneID", "seqname", "start", "end", "gene_info", "logic_name")]

# looks like all the "protein_coding" ones are unhelpful. let's just extract the gene names
nar_annot_results[c("type", "gene")] <- str_split_fixed(nar_annot_results$gene_info, '=',2)

nar_unique <- as.data.frame(nar_annot_results$gene)
colnames(nar_unique) <- "gene"

nar_unique <- subset(nar_unique, gene!="protein_coding")

# Save gene lists
write.table(bow_unique, "bowhead_cand_genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
write.table(nar_unique, "narwhal_cand_genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)

# Also save gene list with annot info (narwhal). bowhead has gene name already in original annot file
write.table(nar_annot_results, "narwhal_cand_genes_annot.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names = TRUE)

### Next, looking at go terms with Enrichr
# Taking a look at what databases are available in EnrichR
listEnrichrDbs()

## Bowhead
# Run enrichR, searching the biological process and molecular function databases
bow_enriched <- enrichr(bow_unique$gene, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

# Pull the results from one of the databases, biological process first
bow_enriched_BP <- as_tibble(bow_enriched[["GO_Biological_Process_2018"]])
bow_enriched_MF <- as_tibble(bow_enriched[["GO_Molecular_Function_2018"]])
bow_enriched_CC <- as_tibble(bow_enriched[["GO_Cellular_Component_2018"]])

## Narwhal
# Run enrichR, searching the biological process and molecular function databases
nar_enriched <- enrichr(nar_unique$gene, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"))

# Pull the results from one of the databases, biological process first
nar_enriched_BP <- as_tibble(nar_enriched[["GO_Biological_Process_2018"]])
nar_enriched_MF <- as_tibble(nar_enriched[["GO_Molecular_Function_2018"]])
nar_enriched_CC <- as_tibble(nar_enriched[["GO_Cellular_Component_2018"]])

# Selecting results with adjusted.p.value < 0.05
bow_enriched_BP_filter <- dplyr::filter(bow_enriched_BP, Adjusted.P.value < 0.05)
bow_enriched_BP_filter$type <- "bowhead_BP"

bow_enriched_MF_filter <- dplyr::filter(bow_enriched_MF, Adjusted.P.value < 0.05)
bow_enriched_MF_filter$type <- "bowhead_MF"

bow_enriched_CC_filter <- dplyr::filter(bow_enriched_CC, Adjusted.P.value < 0.05)
bow_enriched_CC_filter$type <- "bowhead_CC"

nar_enriched_BP_filter <- dplyr::filter(nar_enriched_BP, Adjusted.P.value < 0.05)
nar_enriched_BP_filter$type <- "narwhal_BP"

nar_enriched_MF_filter <- dplyr::filter(nar_enriched_MF, Adjusted.P.value < 0.05)
nar_enriched_MF_filter$type <- "narwhal_MF"

nar_enriched_CC_filter <- dplyr::filter(nar_enriched_CC, Adjusted.P.value < 0.05)
nar_enriched_CC_filter$type <- "narwhal_CC"

all <- rbind(bow_enriched_BP_filter,bow_enriched_MF_filter,bow_enriched_CC_filter,nar_enriched_BP_filter,nar_enriched_MF_filter,nar_enriched_CC_filter)

write.csv(all, "GO_results_adjpvalue0.05_20kb.csv", row.names=F)
