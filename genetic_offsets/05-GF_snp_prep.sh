# Use vcftools to select just the target SNPs for narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.vcf.gz

# (filtering for just the snp list we determined from genome scan and lfmms)

vcftools --gzvcf narwhal_snps.PASS.minQ50.miss.biallel.min100kb.autosomes.n57.maf.miss01.sm.imputed.thin1000.vcf.gz --positions narwhal_top0.01_scan_4env_merged_snplist.txt --recode --recode-INFO-all --stdout > narwhal_top0.01_scan_4env_snps.vcf

# Run vcf2forR.sh to convert to format for GradientForest
./vcf2forR.sh narwhal_top0.01_scan_4env_snps
