#!/bin/bash

# want to make scaffold name match gff file
# basically go from “ENA|SIHG01000001|SIHG01000001.1 Monodon monoceros isolate NGI Contig0_ctg1, whole genome shotgun seq” to just “SIHG01000001.1”

# keep characters 18-31
# want to go by letter but not sure how, workaround is keep characters between the two "|" and then add ".1" to each one to match gff file

awk -F '|' '/^>/ { $0 = ">" $2 } { print }' NAR_GCF_005190385.2.fasta | sed '/^>/s/$/.1/' > NAR_GCF_005190385.2.scafname.fasta