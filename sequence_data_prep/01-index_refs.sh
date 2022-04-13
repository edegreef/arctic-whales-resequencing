#!/bin/bash

# Index references (took about an hour for one 2GB genome)
bwa index ref_genomes/BOW_reference.fasta
samtools faidx ref_genomes/BOW_reference.fasta

bwa index ref_genomes/NAR_reference.fasta
samtools faidx ref_genomes/NAR_reference.fasta
