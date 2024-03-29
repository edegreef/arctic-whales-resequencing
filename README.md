# arctic-whales-resequencing
This is a repository for scripts used in analyzing Canadian Arctic narwhal (*Monodon monoceros*) and bowhead whale (*Balaena mysticetus*) resequencing data. These scripts were executed on computing resources through the University of Manitoba and the Digitial Research Alliance of Canada. 

### Sequencing data preparation [:file_folder:](https://github.com/edegreef/arctic-whales-resequencing/tree/main/sequence_data_prep)
Steps to process raw reads (fastq.gz's) to prepare for SNP calling.
01. Index reference genomes and look at stats
02. Run fastqc on reads for quality check
03. Merge fastqs if multiple lanes. Will note here that it's better to trim first before merging. In this case it didn't make a difference, but in future it's better practice to trim before merge in case of a bad sequencing lane or something.
04. Trim fastqs with *Trimmomatic*
05. Map reads to reference genome with *BWA*. Two scripts here, one for compute canada (05cc) and one for biology cluster (05biol). Did it in parts/batches across both.
06. Remove duplicate reads with *picard*
07. Add read groups with *picard*
08. Remove non-primary alignments with *samtools*
09. Check coverages and count reads
10. Downsample samples with higher coverage with *gatk* (to avoid biases in SNP calling)


### SNP data preparation [:file_folder:](https://github.com/edegreef/arctic-whales-resequencing/tree/main/snps_prep)
Steps to call and filter SNPs
01. Call SNPs with *platypus*
02. Edit sample IDs (remove path from sample IDs and removing sequencing suffix)
03. Look at SNP stats
04. Notes on running DifCover to identify X and Y
05. Look at DifCover results in R
06. Annotate which windows are sex-linked and control using *bedtools*
07. Finalize sex-linked scaffolds
08. SNP filtering pipeline (includes removing indels, QUAL filter, MQ filter, QD filter, missingness, biallelic, removing small scaffolds, removing sex-linked scaffolds, HWE filter, MAF filter, LD pruning...)


### Population structure [:file_folder:](https://github.com/edegreef/arctic-whales-resequencing/tree/main/pop_structure)
* Principcal Component Analysis (PCA) with *pcadapt*: analyzed and plotted with "PCA_pcadapt.R"
* Admixture with sparse Non-Negative Matrix Factorization (sNMF) in *LEA*: analyzed with "SNMF_lea.R", then plotted admixture results in "admixture_plot.R"
* Pairwise differentiation (Reich's Fst): estimated through "FST_Reich.R", then looked at isolation-by-distance with "map_distances_and_IBD.R"
* Runs of Homozygosity (ROH) with *plink*: estimated with "ROH_tadj_pi.sh" then plotted in R with "ROH_tadj_pi_plot.R"
  
### Demographic history [:file_folder:](https://github.com/edegreef/arctic-whales-resequencing/tree/main/demography)
* *PSMC*:
  1. Prep files, including calling SNPs with high coverage bam files
  2. Run PSMC
  3. Run 100 bootstraps
  4. Plot through PSMC
  5. Nicer plot in R
* *SMC++*:
  1. Prep files for SMC++
  2. Run vcf2smc for each scaffold
  3. Run SMC++ estimate with 100 iterations
  4. Plot through SMC++
  5. Nicer plot (combined with PSMC) is in main demography folder
* *GONE*:
  1. Prep files for GONE
  2. Run GONE
  3. Plot GONE in R
* *strataG*
  1. Prep files for strataG
  2. Run ldNe function in strataG
* *EPOS* - extra, exploratory
  1. Prep SFS files
  2. Run epos
  3. Plot results in R

### Adaptation [:file_folder:](https://github.com/edegreef/arctic-whales-resequencing/tree/main/adaptation)
01. Prepare SNPs additional filters, and then impute missing SNPs with *beagle*
02. Prepare haplohh file with *rehh*. Additional script to run R on linux.
03. Look at ihs results and pull out candidate snps with *rehh*. Script includes plot.
04. Prepare input files for next step (SNP and gene locations) 
05. Extract genes within windows of candidate SNPs with *magma*
06. Look at GO-term enrichment with *enrichr*

