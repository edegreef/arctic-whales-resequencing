#!/bin/bash

# 1) First round of filtering, for QUAL, MQ, QD

# Input snp prefix for input and outputs, and reference genome info
snps=BOW_SNPS
path_to_ref=/home/degreefe/whales/ref_genomes/
ref_genome=BOW_reference

# Gatk seemed to run better when going through each filtering criteria one at a time and on command line (crashed if I combined them, not sure why)
# Need reference genome .dict file in same directory as reference fasta
java -jar /home/degreefe/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar CreateSequenceDictionary -R $path_to_ref/$ref_genome.fasta -O $path_to_ref/$ref_genome.dict 

# Filter out QUAL < 50 (QUAL < 30 for Narwhal)
java -jar /home/degreefe/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantFiltration -R $path_to_ref/$ref_genome.fasta -V $snps.vcf.gz -O $snps.QUAL.vcf.gz --filter-name "QUALlt50" --filter-expression "QUAL < 50" 

# Switching to SelectVariants function for next filter steps because the VariantFiltration does not exclude NAs
# Filter out MQ < 40
java -jar /home/degreefe/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -R $path_to_ref/$ref_genome.fasta -V $snps.QUAL.vcf.gz -O $snps.QUAL.MQ.vcf.gz -select "MQ >= 40.0"

# Filter out QD < 4
java -jar /home/degreefe/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -R $path_to_ref/$ref_genome.fasta -V $snps.QUAL.MQ.vcf.gz -O $snps.QUAL.MQ.QD.vcf.gz -select "QD >= 4.0"

# Clean up files and rename final one with "filter1"
mv $snps.QUAL.MQ.QD.vcf.gz $snps.filter1.vcf.gz
mv $snps.QUAL.MQ.QD.vcf.gz.tbi $snps.filter1.vcf.gz.tbi
rm $snps.QUAL*

# 2) Second round of filtering, for max missingness, non-biallelic sites

# Filter out snps with high missingness (in vcftools 1=no missing, 0=all missing)
vcftools --gzvcf $snps.filter1.vcf.gz --max-missing 0.75 --recode --recode-INFO-all --out $snps.filter1.miss

# Remove non-biallelic sites
vcftools --vcf $snps.filter1.miss.recode.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out $snps.filter1.miss.biallel

mv $snps.filter1.miss.biallel.recode.vcf $snps.filter1.miss.biallel.vcf
