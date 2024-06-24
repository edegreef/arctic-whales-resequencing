#!/bin/bash

# Only keep scaffolds minium 10Mb length

# bowhead
vcftools --vcf bowhead_RWmap_snps_rmfix.filter1.miss.biallel.min100kb.autosomes.hwe.n20.vcf --chr CM053039.1 --chr CM053040.1 --chr CM053041.1 --chr CM053042.1 --chr CM053043.1 --chr CM053044.1 --chr CM053045.1 --chr CM053046.1 --chr CM053047.1 --chr CM053048.1 --chr CM053049.1 --chr CM053050.1 --chr CM053051.1 --chr CM053052.1 --chr CM053053.1 --chr CM053054.1 --chr CM053055.1 --chr CM053056.1 --chr CM053057.1 --chr CM053058.1 --recode --recode-INFO-all --out bowhead_RWmap_snps_rmfix.filter1.miss.biallel.min100kb.autosomes.hwe.n20.10MB

# narwhal EBB (BI)
vcftools --vcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.EBB.vcf --chr SIHG01006960.1 --chr SIHG01006964.1 --chr SIHG01006954.1 --chr SIHG01006953.1 --chr SIHG01006963.1 --chr SIHG01006965.1 --chr SIHG01006968.1 --chr SIHG01006958.1 --chr SIHG01006970.1 --chr SIHG01006962.1 --chr SIHG01006966.1 --chr SIHG01006969.1 --chr SIHG01006957.1 --chr SIHG01006961.1 --chr SIHG01006951.1 --chr SIHG01006955.1 --chr SIHG01006967.1 --chr SIHG01006956.1 --chr SIHG01006971.1 --chr SIHG01006952.1 --chr SIHG01006972.1 --chr SIHG01006959.1 --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.EBB.10MB
# narwhal WBB (CHA)
vcftools --vcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.WBB.vcf --chr SIHG01006960.1 --chr SIHG01006964.1 --chr SIHG01006954.1 --chr SIHG01006953.1 --chr SIHG01006963.1 --chr SIHG01006965.1 --chr SIHG01006968.1 --chr SIHG01006958.1 --chr SIHG01006970.1 --chr SIHG01006962.1 --chr SIHG01006966.1 --chr SIHG01006969.1 --chr SIHG01006957.1 --chr SIHG01006961.1 --chr SIHG01006951.1 --chr SIHG01006955.1 --chr SIHG01006967.1 --chr SIHG01006956.1 --chr SIHG01006971.1 --chr SIHG01006952.1 --chr SIHG01006972.1 --chr SIHG01006959.1 --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.WBB.10MB
# narwhal RB (NHB)
vcftools --vcf narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.RB.vcf --chr SIHG01006960.1 --chr SIHG01006964.1 --chr SIHG01006954.1 --chr SIHG01006953.1 --chr SIHG01006963.1 --chr SIHG01006965.1 --chr SIHG01006968.1 --chr SIHG01006958.1 --chr SIHG01006970.1 --chr SIHG01006962.1 --chr SIHG01006966.1 --chr SIHG01006969.1 --chr SIHG01006957.1 --chr SIHG01006961.1 --chr SIHG01006951.1 --chr SIHG01006955.1 --chr SIHG01006967.1 --chr SIHG01006956.1 --chr SIHG01006971.1 --chr SIHG01006952.1 --chr SIHG01006972.1 --chr SIHG01006959.1 --recode --recode-INFO-all --out narwhal_snps.filter1.miss.biallel.min100kb.autosomes.hwe.RB.10MB
