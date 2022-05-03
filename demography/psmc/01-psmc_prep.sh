#!/bin/bash

# install psmc
#git clone https://github.com/lh3/psmc
#cd psmc
#make
#cd utils
#make

# test
#cd ..
#utils/fq2psmcfa -q20 diploid.fq.gz > diploid.psmcfa
#./psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o diploid.psmc diploid.psmcfa
#utils/psmc2history.pl diploid.psmc | utils/history2ms.pl > ms-cmd.sh
#utils/psmc_plot.pl diploid diploid.psmc

#download vcfutils https://github.com/lh3/samtools/blob/master/bcftools/vcfutils.pl
#chmod +x vcfutils.pl


# Bowhead

reference=/home/degreefe/whales/ref_genomes/BOW_reference.fasta
list=scaf_min100kb 
bam=88_Pang_S51.deDupRG.bam
ID=88_Pang

# 88_Pang is 19x

# -d refers to min read depth (1/3 average depth) and -D to max (x2 average depth) .going with avg x9
## this loop takes maybe a day for this batch

while read scaffold
do
samtools mpileup -q 20 -Q 20 -C 50 -u -r $scaffold -f $reference $bam | /home/degreefe/programs/bcftools-1.9/bcftools call -c | ./vcfutils.pl vcf2fq -d 6 -D 38 > prep/$ID.$scaffold.fq
done < $list

# next need to merge the fqs
cd prep
cat $ID*.fq > $ID.consensus.fq

# copy consensus file and go back to previous directory
cp *consensus.fq /home/degreefe/whales/BOW_snps/psmc/
cd ..

# convert fastq to input for PSMC
/home/degreefe/programs/psmc/utils/fq2psmcfa $ID.consensus.fq > $ID.psmcfa