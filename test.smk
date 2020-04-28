# snakejob -s test.smk -j 100
# snakemake -s test.smk
import os.path
from os import path
import pandas as pd

shell.prefix('module load samtools singularity freebayes/1.3.2 vcflib bcftools tabix; ')

configfile: "test.yaml"
SAMPLES = pd.read_csv("test.csv", index_col='SampleID')

rule all:
    input:
        expand("data/test/{SAMPLE_IDS}.chrM.bam", SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("data/test/{SAMPLE_IDS}.chrM.bam.bai", SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("data/test/{SAMPLE_IDS}MTData.txt", SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("data/test/{SAMPLE_IDS}.chrM.fb.vcf.gz", SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("data/test/{SAMPLE_IDS}.haplogrep.txt", SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("data/test/{SAMPLE_IDS}.phymer.txt", SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("data/test/{SAMPLE_IDS}.chrM.mity.vcf.gz", SAMPLE_IDS = SAMPLES.index.tolist())


rule copy:
    input:
        bam = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['wgs'],
        bai = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['index']
    output:
        bam = 'data/test/{SAMPLE_IDS}.chrM.bam',
        bai = 'data/test/{SAMPLE_IDS}.chrM.bam.bai'
    run:
        ## For samples aligned to hg37 extract chromsome MT
        if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37:
            shell('samtools view -b {input.bam} MT > {output.bam}; samtools index -b {output.bam}')
            print("Using hg37 as reference build")
        ## For samples aligned to hg38 extract chromsome chrM
        elif SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38:
            shell('samtools view -b {input.bam} chrM > {output.bam}; samtools index -b {output.bam}')
            print("Using hg38 as reference build")
        else:
            print(SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] + ': Build Not Recognised')

rule fastMitoCalc:
    input:
        bam = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['wgs'],
        bai = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['index']
    output: "data/test/{SAMPLE_IDS}MTData.txt"
    params:
        wd = 'data/test'
    run:
        ## For samples aligned to hg37 use default chromosome prefix (none) and mtname (MT)
        if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37:
            shell('perl src/fastMitoCalc/fastMitoCalc.pl -f {input.bam} -w {params.wd} -p src/fastMitoCalc/BaseCoverage')
            print("Using hg37 as reference build")
        ## For samples aligned to hg38 specify chromosome prefix (chr) and mtname (chrM)
        elif SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38:
            shell('perl src/fastMitoCalc/fastMitoCalc.pl -f {input.bam} -e chr -m chrM -w {params.wd} -p src/fastMitoCalc/BaseCoverage')
            print("Using hg38 as reference build")
        else:
            print(SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] + ': Build Not Recognised')

rule mity:
    input:
        bam = 'temp/{SAMPLE_IDS}.chrM.bam',
        bai = 'temp/{SAMPLE_IDS}.chrM.bam.bai'
    output:
        vcf = "data/test/{SAMPLE_IDS}.chrM.mity.vcf.gz",
        tbi = "data/test/{SAMPLE_IDS}.chrM.mity.vcf.gz.tbi",
    params:
        ref = lambda wildcards: 'hs37d5' if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37 else ('hg38' if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38 else None),
        wd = 'data/test'
    shell:
        'singularity run src/mity.sif call --reference {params.ref} --out-folder-path {params.wd} --normalise {input.bam}'


rule freebayes:
    input:
        bam = 'temp/{SAMPLE_IDS}.chrM.bam',
        bai = 'temp/{SAMPLE_IDS}.chrM.bam.bai',
        fa = lambda wildcards: 'raw/reference/hs37d5.MT.fa' if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37 else ('raw/reference/hg38.chrM.fa' if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38 else None)
    output:
        vcf = "data/test/{SAMPLE_IDS}.chrM.fb.vcf",
    shell:
        """
        freebayes -f {input.fa} -b {input.bam} \
        --min-mapping-quality 30 \
        --min-base-quality 24 \
        --min-alternate-fraction 0.6 \
        --min-alternate-count 4 \
        --ploidy 1 \
        | vcffilter -f "QUAL > 20" > {output.vcf}
        """

rule bgzip:
    input: vcf = "data/test/{SAMPLE_IDS}.chrM.fb.vcf"
    output:
        bcf = "data/test/{SAMPLE_IDS}.chrM.fb.vcf.gz",
        tbi = "data/test/{SAMPLE_IDS}.chrM.fb.vcf.gz.tbi"
    shell:
        'bgzip {input.vcf}; tabix {output.bcf}'


rule haplogrep:
    input:
#        vcf = "data/mity/{SAMPLE_IDS}.chrM.mity.vcf.gz",
        vcf = "data/test/{SAMPLE_IDS}.chrM.fb.vcf.gz",
    output:
        out = "data/test/{SAMPLE_IDS}.haplogrep.txt",
    shell:
        'java -jar src/haplogrep-2.1.25.jar --in {input.vcf} --format vcf --out {output.out}'


rule phymer:
    input:
        bam = 'temp/{SAMPLE_IDS}.chrM.bam',
        bai = 'temp/{SAMPLE_IDS}.chrM.bam.bai'
    output:
        out = "data/test/{SAMPLE_IDS}.phymer.txt",
    shell:
        """
        ./src/phymer/Phy-Mer.py --print-ranking --verbose \
        --min-DoC=10 src/phymer/PhyloTree_b16_k12.txt {input.bam} > {output.out}
        """
