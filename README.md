# arctic-whales-resequencing
Currently a working repository


sequence_data_prep - steps from raw sequence fastq.gz files to snp vcf

- remember that it's better to trim before merging fastq.gz! in the case of narwhal and bowhead dataset it didn't make difference, but can make difference for other situations esp if 1 lane is bad-- better practice to to trim then merge

snps_prep - steps for snps filtering (including identifying sex-linked scaffolds to filter)

demography - running smc++, epos, & psmc
