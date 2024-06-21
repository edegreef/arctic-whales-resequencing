# remove sites that are same AA, RR, or missmiss across all samples

/home/degreefe/programs/bcftools-1.9/bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES || COUNT(GT="mis")=N_SAMPLES'  bowhead_RWmap_allvariants.ID.vcf.gz > bowhead_RWmap_allvariants.ID.removefix1.vcf

bgzip bowhead_RWmap_allvariants.ID.removefix1.vcf

# Can also run fixedsnps.R to identify ones that also match but have some missing data, but these end up getting filtered anyway in the downstream filters.
# If using fixedsnps.R then use "fixed_snps_CHROM_POS_part2_allvariants.txt" to remove the last few fixed snps with missing data
#vcftools --gzvcf bowhead_RWmap_allvariants.ID.removefix1.vcf.gz --exclude-positions fixed_snps_CHROM_POS_part2_allvariants.txt --recode --recode-INFO-all --out bowhead_RWmap_allvariants.ID.rmfix
