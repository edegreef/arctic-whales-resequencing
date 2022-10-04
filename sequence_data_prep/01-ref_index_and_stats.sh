#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --job-name=reference_stuff
#SBATCH --mem=10GB
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@myucdavis.edu
#SBATCH --output=%x-%j.out

# load modules
module load StdEnv/2020 bwa/0.7.17 samtools/1.15.1

#1) Index references (took about an hour for one 2GB genome)
# bowhead
bwa index ref_genomes/BOW_reference.fasta
samtools faidx ref_genomes/BOW_reference.fasta

# narwhal
bwa index ref_genomes/NAR_GCF_005190385.2.scafname.fasta
samtools faidx ref_genomes/NAR_GCF_005190385.2.scafname.fasta

#2) Assemblathon (need FAlite.pm in working directory, downloaded from KorfLab github page along with assemblathon_stats.pl script -- https://github.com/ucdavis-bioinformatics/assemblathon2-analysis)
# this used to work easily, seems like may need older perl version
# also to bypass the qw(...) error in assemblathon_stats.pl, can manually add parentheses for lines 301 & 428 (qw(...))
module load nixpkgs/16.09 perl/5.16.3
perl assemblathon_stats.pl ref_genomes/BOW_reference.fasta > BOW_reference.assemblathon.txt
perl assemblathon_stats.pl ref_genomes/NAR_GCF_005190385.2.scafname.fasta > NAR_GCF_005190385.2.scafname.assemblathon.txt

#3) Extract scaffold lengths
cat ref_genomes/BOW_reference.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > BOW_reference.scaffoldlengths.csv
cat ref_genomes/NAR_GCF_005190385.2.scafname.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > NAR_GCF_005190385.2.scafname.scaffoldlengths.csv

