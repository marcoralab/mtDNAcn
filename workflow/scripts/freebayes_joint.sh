#!/bin/sh
bsub -J "FB.ridge" -P acc_load -q premium -n 24 -R span[hosts=1] -R himem -R rusage[mem=50000] -W 144:00 -L /bin/bash \
-o data/freebayes/joint/FB.ridge.stdout -eo data/freebayes/joint/FB.ridge.stderr \
"module load freebayes/1.3.2 vcflib; freebayes -f raw/reference/hs37d5.MT.fa -L MAP_chrMT_bamFiles.txt --min-mapping-quality 30 --min-base-quality 24 --min-alternate-fraction 0.6 --min-alternate-count 4 --use-best-n-alleles 3 --ploidy 1 | vcffilter -f \"QUAL > 20\" > data/freebayes/joint/MAP.chrmMT.joint.vcf"
