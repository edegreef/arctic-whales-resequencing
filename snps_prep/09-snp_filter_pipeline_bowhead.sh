#!/bin/bash

# Script for doing the rest of the SNP filters.
# Including:
# Step 1: Indel and PASS filter
# Step 2: QUAL filter, MQ filter, QD filter (after this step, renamed to "fitler1")
# Step 3: Missingness filter, biallelic filter
# Step 4: Removing small scaffolds (<100kb)
# Step 5: Removing sex-linked snps
# Step 6: Additional HWE, MAF filter
# Step 7: LD pruning


# Make sure edit ID is done to vcf beforehand
# Also need filter.hets.R file in directory for the HWE filter to work

# Set up these before running
variant_prefix=bowhead_RWmap_allvariants.ID.rmfix
snps_prefix=bowhead_RWmap_snps_rmfix
bcftools_path=/home/degreefe/programs/bcftools-1.9
gatk_jar_path=/home/degreefe/programs/gatk-4.1.9.0
ref_genome_path=/home/degreefe/whales/ref_genomes
ref_genome=GCA_028564815.2_mEubGla1.1.hap2._XY_genomic.fna
min_scaf_list=scaf_min100kb_RW #list of scaffolds to keep
plink_path=/home/degreefe/programs/



######## Step 1: Remove indels and SNPs that don't have the "PASS" filter. Start with vcf.gz file (with corresponding tbi)
# Remove indels
vcftools --gzvcf $variant_prefix.vcf.gz --remove-indels --recode --recode-INFO-all --out $snps_prefix

mv $snps_prefix.recode.vcf $snps_prefix.vcf

# Zip & index
bgzip $snps_prefix.vcf
tabix -p vcf $snps_prefix.vcf.gz

# Filter out non-PASS sites
$bcftools_path/bcftools view -f PASS $snps_prefix.vcf.gz -o $snps_prefix.PASS.vcf

# Zip & index
bgzip $snps_prefix.PASS.vcf
tabix -p vcf $snps_prefix.PASS.vcf.gz



######## Step 2: Next filtering round
# Need reference genome .dict file in same directory as reference fasta
#java -jar /home/degreefe/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar CreateSequenceDictionary -R $ref_genome_path/$ref_genome -O $ref_genome_path/GCA_028564815.2_mEubGla1.1.hap2._XY_genomic.dict 

# Filter out QUAL < 50
java -jar $gatk_jar_path/gatk-package-4.1.9.0-local.jar VariantFiltration -R $ref_genome_path/$ref_genome -V $snps_prefix.PASS.vcf.gz -O $snps_prefix.PASS.QUAL.vcf.gz --filter-name "QUALlt50" --filter-expression "QUAL < 50" 

# Filter out MQ < 40
java -jar $gatk_jar_path/gatk-package-4.1.9.0-local.jar SelectVariants -R $ref_genome_path/$ref_genome -V $snps_prefix.PASS.QUAL.vcf.gz -O $snps_prefix.PASS.QUAL.MQ.vcf.gz -select "MQ >= 40.0"

# Filter out QD < 4
java -jar $gatk_jar_path/gatk-package-4.1.9.0-local.jar SelectVariants -R $ref_genome_path/$ref_genome -V $snps_prefix.PASS.QUAL.MQ.vcf.gz -O $snps_prefix.PASS.QUAL.MQ.QD.vcf.gz -select "QD >= 4.0"

# Rename with "filter1" and clean up files
mv $snps_prefix.PASS.QUAL.MQ.QD.vcf.gz $snps_prefix.filter1.vcf.gz
mv $snps_prefix.PASS.QUAL.MQ.QD.vcf.gz.tbi $snps_prefix.filter1.vcf.gz.tbi

rm $snps_prefix.PASS.QUAL.vcf.gz $snps_prefix.PASS.QUAL.vcf.gz.tbi $snps_prefix.PASS.QUAL.MQ.vcf.gz $snps_prefix.PASS.QUAL.MQ.vcf.gz.tbi


######## Step 3: Missingness & bi-allelic filter. 

# Filter out snps with high missingness (in vcftools 1=no missing, 0=all missing) and remove non-biallelic snps
vcftools --gzvcf $snps_prefix.filter1.vcf.gz --max-missing 0.75 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out $snps_prefix.filter1.miss.biallel

mv $snps_prefix.filter1.miss.biallel.recode.vcf $snps_prefix.filter1.miss.biallel.vcf

# Zip & index
bgzip $snps_prefix.filter1.miss.biallel.vcf
tabix -p vcf $snps_prefix.filter1.miss.biallel.vcf.gz

######## Step 4: Remove scaffolds less than 100kb in length. Needs input list of scaffold names to KEEP. I chose 100kb based on my dataset, but should be checked/adjusted for other data (example, if this filter loses a very large chunk of the genome, may need to use a smaller scaffold filter like 10kb).

# Convert scaffold list to one-liner to use in bcftools
awk '{print $1}' $min_scaf_list | paste -s -d, - > scaf_list_line

# Set up list
list=`cat scaf_list_line`

# Filter vcf for these scaffolds
$bcftools_path/bcftools filter --regions $list $snps_prefix.filter1.miss.biallel.vcf.gz > $snps_prefix.filter1.miss.biallel.min100kb.vcf

# Can make a list of scaffolds to double check
grep -v "^#" $snps_prefix.filter1.miss.biallel.min100kb.vcf | cut -f1 | sort | uniq > filtered_contig_list_check.txt

# Zip & index
bgzip $snps_prefix.filter1.miss.biallel.min100kb.vcf
tabix -p vcf $snps_prefix.filter1.miss.biallel.min100kb.vcf.gz

######## Step 5: Filter out sex-linked SNPS to create autosomal dataset
# easy for right whale mapped because have X and Y:
vcftools --gzvcf $snps_prefix.filter1.miss.biallel.min100kb.vcf.gz --not-chr CM053059.1 --not-chr CM053080.1 --recode --recode-INFO-all --out $snps_prefix.filter1.miss.biallel.min100kb.autosomes

mv $snps_prefix.filter1.miss.biallel.min100kb.autosomes.recode.vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.vcf

bgzip $snps_prefix.filter1.miss.biallel.min100kb.autosomes.vcf
tabix -p vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.vcf.gz



######## Step 6: Additional HWE, MAF filter
# Remove SNPs out of Hardy-Weinberg Equilibrium
# First need to create the out.hwe file
vcftools --gzvcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.vcf.gz --hardy

# Run filter.hets.R script to create list of snps with het > 0.6
R --vanilla < filter.hets.R
vcftools --gzvcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.vcf.gz --exclude-positions snps_het06_filter --recode --recode-INFO-all --out $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe

mv $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.recode.vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.vcf
bgzip $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.vcf
tabix -p vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz

# Minor-allele frequency filter
vcftools --gzvcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.vcf.gz --maf 0.05 --recode --recode-INFO-all --out $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf

mv $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.recode.vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf

bgzip $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf
tabix -p vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf.gz

######## Step 7: LD pruning

# Create list of snps to prune
# --double-id b/c multiple "_" in sample name
$plink_path/plink --vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf.gz --allow-extra-chr --indep-pairwise 50 5 0.8 --set-missing-var-ids @:#\$1,\$2 --out pruned.50window.r0.8 --double-id

# Preparing the prune-in file to filter vcf
sed 's/:/\t/g' pruned.50window.r0.8.prune.in > temp.in
sed 's/...$//' temp.in > temp2.in
mv temp2.in list.pruned.50window.r0.8.prune.in
rm temp.in

# Prune the snps
vcftools --gzvcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.vcf.gz --positions list.pruned.50window.r0.8.prune.in --recode --recode-INFO-all --out $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08

mv $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.recode.vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.vcf

bgzip $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.vcf
tabix -p vcf $snps_prefix.filter1.miss.biallel.min100kb.autosomes.hwe.maf.LDprunedr08.vcf.gz

