# Script to prepare file formating for magma, using genome annotations for narwhal and bowhead

# GENELOC = gene locations (4 columns: gene ID, chrom, start site, stop site)

library(tidyverse)
setwd("C:/Users/eveli/Dropbox/Whales with Garroway/ref_genomes/annotations/")

# Start from gff file

### BOWHEAD WHALE:
# Load annotation
gff <- read.table("BOW_genomegff/B_mysticetus_cleaned.gff") #bowhead

# Add header/col names (https://useast.ensembl.org/info/website/upload/gff.html)
colnames(gff) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Need to extract ID from "attribute" column
gff[c("attribute_ID", "gene_name")] <- str_split_fixed(gff$attribute, ';',2)

# Later realizing magma doesn't like duplicates, so going to add the hit number too 
gff[c("scaffold", "hit", "hit_number", "0_0")] <- str_split_fixed(gff$attribute_ID, ':',4)

# Add scaffold name too for the unique ID because a few "duplicates" but on different scaffolds
gff$gene_ID <- paste(gff$seqname, gff$hit_number, gff$gene_name, sep="-")

# Prepare the gene loc file
geneloc <- gff[c("gene_ID", "seqname", "start", "end")]

# Check for duplicates
dups <- geneloc[duplicated(geneloc$gene_ID)|duplicated(geneloc$gene_ID, fromLAST=TRUE),]

# Save file
write.table(geneloc, "bowhead_geneloc.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)



### NARWHAL:
# Load annotation
# Narwhal gff has lots of rows with ## so reading it a different way
gff <- read_delim("Monodon_monoceros.NGI_Narwhal_1.107.gff3", delim="\t",comment="#", col_names = F)

# Add header/col names (https://useast.ensembl.org/info/website/upload/gff.html)
colnames(gff) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Narwhal annotation needs to filter a bit:
features <- unique(gff$feature)
features
gff_genes_only <- subset(gff, feature=="gene")

# Narwhal gff actually has a gene ID so let's use that
gff_genes_only[c("gene_ID", "gene_info", "gene_ID2", "logic_name", "version")] <- str_split_fixed(gff_genes_only$attribute, ';',5)

# Just need the gene IDs starting with ENSMMNG in the first split column
# Split "gene_ID"
gff_genes_only[c("ID", "ENSMMNG")] <- str_split_fixed(gff_genes_only$gene_ID, ':',2)

# Prepare the gene loc file
geneloc <- gff_genes_only[c("ENSMMNG", "seqname", "start", "end")]

# Save file
write.table(geneloc, "narwhal_geneloc.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)

