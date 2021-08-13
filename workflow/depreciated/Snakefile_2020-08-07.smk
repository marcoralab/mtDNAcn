## Snakefile for calculating mtDNA copy number

import os.path
from os import path
import pandas as pd

shell.prefix('module load samtools singularity freebayes/1.3.2 vcflib bcftools tabix; ')

configfile: "config/config.yaml"
SAMPLES = pd.read_csv("sandbox/SampleInfo_test.csv", index_col='SampleID')
#SAMPLES = pd.read_csv("data/SampleInfo.csv", index_col='SampleID')

rule all:
    input:
        expand('temp/{SAMPLE_IDS}.chrM.bam', SAMPLE_IDS = SAMPLES.index.tolist()),
        # expand("data/mitocalc/{SAMPLE_IDS}MTData.txt", SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("data/mity/{SAMPLE_IDS}.chrM.mity.vcf.gz", SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("data/phymer/{SAMPLE_IDS}.phymer.txt", SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("data/haplogrep/{SAMPLE_IDS}.haplogrep.txt", SAMPLE_IDS  = SAMPLES.index.tolist()),
        expand("data/freebayes/{SAMPLE_IDS}.chrM.fb.vcf.gz", SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("data/mosdepth/{SAMPLE_IDS}.mosdepth.summary.txt", SAMPLE_IDS = SAMPLES.index.tolist()),
        # 'data/freebayes/joint/AMPAD.chrM.fb_joint.vcf'

rule copy:
    input:
        bam = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['wgs'],
        bai = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['index']
    output:
        bam = temp('temp/{SAMPLE_IDS}.chrM.bam'),
        bai = temp('temp/{SAMPLE_IDS}.chrM.bam.bai')
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

# rule fastMitoCalc:
#     input:
#         bam = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['wgs'],
#         bai = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['index']
#     output: "data/mitocalc/{SAMPLE_IDS}MTData.txt"
#     params:
#         wd = 'data/mitocalc'
#     run:
#         ## For samples aligned to hg37 use default chromosome prefix (none) and mtname (MT)
#         if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37:
#             shell('perl src/fastMitoCalc/fastMitoCalc.pl -f {input.bam} -w {params.wd} -p src/fastMitoCalc/BaseCoverage')
#             print("Using hg37 as reference build")
#         ## For samples aligned to hg38 specify chromosome prefix (chr) and mtname (chrM)
#         elif SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38:
#             shell('perl src/fastMitoCalc/fastMitoCalc.pl -f {input.bam} -e chr -m chrM -w {params.wd} -p src/fastMitoCalc/BaseCoverage')
#             print("Using hg38 as reference build")
#         else:
#             print(SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] + ': Build Not Recognised')

rule mosdepth:
    input:
        bam = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['wgs'],
        bai = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['index'],
        fa = 'raw/reference/Homo_sapiens_assembly19.fasta' # lambda wildcards: 'raw/reference/Homo_sapiens_assembly19.fasta' if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37 else ('raw/reference/Homo_sapiens_assembly38.fasta' if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38 else None)
    output: "data/mosdepth/{SAMPLE_IDS}.mosdepth.summary.txt"
    params:
        wd = 'data/mosdepth/{SAMPLE_IDS}'
    run:
        ## For samples in bam format
        if SAMPLES.loc[wildcards.SAMPLE_IDS]['ext'] == 'bam':
            shell('./src/mosdepth -n -t 4 {params.wd} {input.bam}')
            print("Calculating coverage from .bam file")
        ## For samples in cram format
        elif SAMPLES.loc[wildcards.SAMPLE_IDS]['ext'] == 'cram':
            shell('./src/mosdepth -n  -t 4 -f {input.fa} {params.wd} {input.bam}')
            print("Calculating coverage from .cram file")
        else:
            print(SAMPLES.loc[wildcards.SAMPLE_IDS]['ext'] + ': File Format not recongised')

rule mity:
    input:
        bam = 'temp/{SAMPLE_IDS}.chrM.bam',
        bai = 'temp/{SAMPLE_IDS}.chrM.bam.bai'
    output:
        vcf = "data/mity/{SAMPLE_IDS}.chrM.mity.vcf.gz",
        tbi = "data/mity/{SAMPLE_IDS}.chrM.mity.vcf.gz.tbi",
    params:
        ref = lambda wildcards: 'hs37d5' if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37 else ('hg38' if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38 else None),
        wd = 'data/mity'
    shell:
        'singularity run src/mity.sif call --reference {params.ref} --out-folder-path {params.wd} --normalise {input.bam}'

rule freebayes:
    input:
        bam = 'temp/{SAMPLE_IDS}.chrM.bam',
        bai = 'temp/{SAMPLE_IDS}.chrM.bam.bai',
        fa = lambda wildcards: 'raw/reference/hs37d5.MT.fa' if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37 else ('raw/reference/hg38.chrM.fa' if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38 else None)
    output:
        vcf = "data/freebayes/{SAMPLE_IDS}.chrM.fb.vcf",
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
    input: vcf = "data/freebayes/{SAMPLE_IDS}.chrM.fb.vcf"
    output:
        bcf = "data/freebayes/{SAMPLE_IDS}.chrM.fb.vcf.gz",
        tbi = "data/freebayes/{SAMPLE_IDS}.chrM.fb.vcf.gz.tbi"
    shell:
        'bgzip {input.vcf}; tabix {output.bcf}'


def bcf_input(wildcards):
    return expand("data/freebayes/{SAMPLE_IDS}.chrM.fb.vcf.gz", SAMPLE_IDS = SAMPLE_IDS)
#
#
# rule vcf_merge:
#     input: vcfs = bcf_input
#     output: vcf = 'data/freebayes/joint/AMPAD.chrM.fb_joint.vcf'
#     shell: '{input.vcfs} > data/freebayes/joint/chrMT_vcfiles.txt'
#         #
#         # """
#         # bcftools merge {input.vcfs} --missing-to-ref | \
#         # bcftools view --types snps | \
#         # bcftools view --exclude-types mnps --max-alleles 2 -Ov -o {output.vcf}'
#         # """

rule haplogrep:
    input:
#        vcf = "data/mity/{SAMPLE_IDS}.chrM.mity.vcf.gz",
        vcf = "data/freebayes/{SAMPLE_IDS}.chrM.fb.vcf.gz",
    output:
        out = "data/haplogrep/{SAMPLE_IDS}.haplogrep.txt",
    shell:
        'java -jar src/haplogrep-2.1.25.jar --in {input.vcf} --format vcf --out {output.out}'

rule phymer:
    input:
        bam = 'temp/{SAMPLE_IDS}.chrM.bam',
        bai = 'temp/{SAMPLE_IDS}.chrM.bam.bai'
    output:
        out = "data/phymer/{SAMPLE_IDS}.phymer.txt",
    shell:
        """
        ./src/phymer/Phy-Mer.py --print-ranking --verbose \
        --min-DoC=10 src/phymer/PhyloTree_b16_k12.txt {input.bam} > {output.out}
        """
