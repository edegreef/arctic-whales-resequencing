#!/bin/bash

# Merge raw fastq.gz files for reads across multiple lanes.
# Code obtained from Paul at https://www.biostars.org/p/317385/

# will add note that if doing on trimmed reads that have "paired" in the file name, might need to adjust a bit
# for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 25- | rev; done | sort | uniq)

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)

    do echo "Merging R1"

cat "$i"_L00*_R1_001.fastq.gz > "$i"_R1_001.fastq.gz

       echo "Merging R2"

cat "$i"_L00*_R2_001.fastq.gz > "$i"_R2_001.fastq.gz

done
