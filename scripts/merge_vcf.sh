#!/bin/sh
bcftools merge data/freebayes/*.vcf.gz --missing-to-ref | \
bcftools view --types snps | \
bcftools view --exclude-types mnps --max-alleles 2 -Ov -o data/freebayes/joint/AMPAD.chrM.fb_joint.vcf; \
java -jar src/haplogrep-2.1.25.jar --in data/freebayes/joint/AMPAD.chrM.fb_joint.vcf --format vcf --out data/freebayes/joint/AMPAD.chrM.haplogrep.txt
