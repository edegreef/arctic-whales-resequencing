![Logo](https://github.com/edegreef/arctic-whales-resequencing/assets/49288304/21b4781f-a1cc-4a05-8703-51f74cafbdb3)

This is a repository for scripts used in analyzing Canadian Arctic narwhal (*Monodon monoceros*) and bowhead whale (*Balaena mysticetus*) resequencing data. These scripts were executed on computing resources through the University of Manitoba and the Digitial Research Alliance of Canada. 

### Sequencing data preparation [:file_folder:](https://github.com/edegreef/arctic-whales-resequencing/tree/main/sequence_data_prep)
Steps to process raw reads (fastq.gz's) to prepare for SNP calling.
01. Index reference genomes and look at stats
02. Run fastqc on reads for quality check
03. Merge fastqs if multiple lanes. Will note here that it's better to trim first before merging. In this case it didn't make a difference, but in future it's better practice to trim before merge in case of a bad sequencing lane.
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
08. Removing fixed SNPs in bowhead whale (because mapped to a different species)
09. SNP filtering pipeline (includes removing indels, QUAL filter, MQ filter, QD filter, missingness, biallelic, removing small scaffolds, removing sex-linked scaffolds, HWE filter, MAF filter, LD pruning...). Separate scripts here for bowhead and narwhal because of slight difference in removing X and Y. The right whale genome used for the bowhead had X and Y chromosomes already identified.
10. A few extra SNP conversions and filter to prepare files into bfile format (bed/bim/fam), check kinships, and remove necessary samples, and then converting to map/ped format.


### Population structure [:file_folder:](https://github.com/edegreef/arctic-whales-resequencing/tree/main/pop_structure)
* Principcal Component Analysis (PCA) with *pcadapt*: analyzed and plotted with "PCA_pcadapt.R"
* Admixture with sparse Non-Negative Matrix Factorization (sNMF) in *LEA*: analyzed with "SNMF_lea.R", then plotted admixture results in "admixture_plot.R"
* Pairwise differentiation (Reich's Fst): estimated through "FST_Reich.R", then looked at isolation-by-distance with "map_distances_and_IBD.R"
* Genome-wide proportion of observed heterozygossity estimated with *ANGSD* to include both multivariant and monomorphic sites. "angsd_heterozygosity_narwhal.sh" and "angsd_heterozygosity_bowhead.sh"
* Proportion of observed heterozygosity estimated with "calculate_het.R" using output from vcftools. This method was done on SNPs and does not include monomorphic sites.
* Runs of Homozygosity (ROH) with *plink*: estimated with "estimate_ROH.sh" then plotted in R with "plot_ROH.R"
  
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
  3. Run SMC++ estimate with 100 iterations/arrays
  4. Plot through SMC++
  5. Nicer plot (combined with PSMC) is in main demography folder "smcpp_and_psmc_plot.R"
  6. Additional scripts for looking at population splits (twopop)
* *GONE*:
  1. Keep only scaffolds >10Mb
  2. Prep files for GONE
  3. Run GONE
  4. Plot GONE in R
* *strataG*
  1. Prep files for strataG
  2. Run ldNe function in strataG
* *EPOS* - extra, exploratory
  1. Prep SFS files
  2. Run epos
  3. Plot results in R
  
