#!/bin/bash

for i in *.bam
do
echo "Indexing: "$i        
samtools index $i
done
