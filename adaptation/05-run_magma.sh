# Magma doesn't like non-numerical scaffold names, so doing a lil editing for format on input files (snploc, geneloc)

# For bowhead files want to remove "scaffold_":
awk '{ gsub(/scaffold_/,"", $2); print } ' bowhead_ihs_cands_snps_nextrmark3_snploc.txt > bowhead_ihs_cands_snps_nextrmark3_snploc_rename.txt
awk '{ gsub(/scaffold_/,"", $2); print } ' bowhead_geneloc.txt > bowhead_geneloc_rename.txt

# Note that magma will not read "0" as a chr (regarding "scaffold_0"). There are no candidate positions on scaffold_0 so it's OK to ignore this here b/c not one of our target regions, but if want to look any genes things later on scaffold_0 will need to keep this in mind.

# For narwhal files want to remove "SIHG0" and the ".1" at the end (not very efficient below, but it works):
awk '{gsub(/SIHG0/,"", $2); print } ' narwhal_ihs_cands_snps_nextrmark3_snploc.txt > narwhal_snploc_temp
awk '{print substr($2,1, (length($2)-2))}' narwhal_snploc_temp > narwhal_cands_col2
awk '{print $1}' narwhal_ihs_cands_snps_nextrmark3_snploc.txt > narwhal_cands_col1
awk '{print $3}' narwhal_ihs_cands_snps_nextrmark3_snploc.txt > narwhal_cands_col3
paste narwhal_cands_col1 narwhal_cands_col2 narwhal_cands_col3 > narwhal_ihs_cands_snps_nextrmark3_snploc_rename.txt

awk '{gsub(/SIHG0/,"", $2); print } ' narwhal_geneloc.txt > narwhal_geneloc_temp
awk '{print substr($2,1, (length($2)-2))}' narwhal_geneloc_temp > narwhal_gene_col2
awk '{print $1}' narwhal_geneloc.txt > narwhal_gene_col1
awk '{print $3, $4}' narwhal_geneloc.txt > narwhal_gene_col34
paste narwhal_gene_col1 narwhal_gene_col2 narwhal_gene_col34 > narwhal_geneloc_rename.txt

# Run magma
# Bowhead whale
/home/degreefe/programs/magma --annotate nonhuman window=20,20 --snp-loc bowhead_ihs_cands_snps_nextrmark3_snploc_rename.txt  --gene-loc bowhead_geneloc_rename.txt --out bowhead_magma_annotation

# Narwhal
/home/degreefe/programs/magma --annotate nonhuman window=20,20 --snp-loc narwhal_ihs_cands_snps_nextrmark3_snploc_rename.txt  --gene-loc narwhal_geneloc_rename.txt --out narwhal_magma_annotation
