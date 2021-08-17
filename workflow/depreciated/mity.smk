import os
from itertools import product
import pandas as pd
RWD = os.getcwd()

shell.prefix('module load singularity; ')
SAMPLES = pd.read_csv('data/SampleInfo.csv', index_col='SampleID')

# ID = ["SM-CJGNJ", "SM-CJEIA", "71992"]
# ID = ["SM-CJGNJ"]

# snakejob -s mity.smk -j 1000 --use-conda --max-jobs-per-second 2
rule all:
    input:
        'data/mity/results/contamination/contamination.txt',
        'data/mity/results/haplogroups/haplogroups.txt',
        'data/mity/results/report/report.html',
        'data/mity/results/haplocheck_res.txt',
        'data/mity/results/haplocheck_res.raw.txt',
        'data/mity/results/haplocheck_res.html'

# rule mity:
#     input:
#         bam = 'temp/{SAMPLE_IDS}.chrM.bam',
#         bai = 'temp/{SAMPLE_IDS}.chrM.bam.bai'
#     output:
#         vcf = "sandbox/mity/{SAMPLE_IDS}.chrM.mity.vcf.gz",
#         tbi = "sandbox/mity/{SAMPLE_IDS}.chrM.mity.vcf.gz.tbi",
#     params:
#         # ref = lambda wildcards: 'hs37d5' if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37 else ('hg38' if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38 else None),
#         ref = "hs37d5",
#         wd = 'sandbox/mity'
#     # envmodules:
#     #     '/hpc/packages/minerva-centos7/singularity/3.2.1/bin/singularity'
#     container:
#         # "docker://joseespinosa/docker-r-ggplot2"
#         'docker://wzhou88/saige:0.44'
#         # "docker://drmjc/mity:0.1.3"
#     shell:
#         'mity call --reference {params.ref} --out-folder-path {params.wd} --region MT:1-500 --normalise {input.bam}'

rule mity:
    input:
        bam = 'temp/{SAMPLE_IDS}.chrM.bam',
        bai = 'temp/{SAMPLE_IDS}.chrM.bam.bai'
    output:
        vcf = "data/mity/calls/{SAMPLE_IDS}.chrM.mity.vcf.gz",
        tbi = "data/mity/calls/{SAMPLE_IDS}.chrM.mity.vcf.gz.tbi",
    params:
        # ref = lambda wildcards: 'hs37d5' if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37 else ('hg38' if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38 else None),
        ref = "hs37d5",
        wd = 'data/mity/calls/'
    # conda:
    #     'workflow/envs/singularity.yaml'
    shell:
        'singularity run src/mity-0.1.3.simg call --reference {params.ref} --out-folder-path {params.wd} --normalise {input.bam}'

def merge_vcf_input(wildcards):
    return expand("data/mity/calls/{SAMPLE_IDS}.chrM.mity.vcf.gz", SAMPLE_IDS = SAMPLES.index.tolist()) # SAMPLES.index.tolist()
def merge_index_input(wildcards):
    return expand("data/mity/calls/{SAMPLE_IDS}.chrM.mity.vcf.gz.tbi", SAMPLE_IDS = SAMPLES.index.tolist()) # SAMPLES.index.tolist()

rule merge_vcfs:
    input:
        vcf = merge_vcf_input,
        index = merge_index_input
    output: "data/mity/results/mity_all.vcf.gz"
    conda:
        'workflow/envs/bcftools.yaml'
    shell: "bcftools merge --missing-to-ref -Oz -o {output} {input.vcf}"

rule haplocheck_hg:
    input: "data/mity/results/mity_all.vcf.gz"
    output:
        'data/mity/results/contamination/contamination.txt',
        'data/mity/results/haplogroups/haplogroups.txt',
        'data/mity/results/report/report.html'
    params:
        out = RWD+'/data/mity/results',
        input = RWD+"/data/mity/results/mity_all.vcf.gz"
    shell: "./src/haplocheck/cloudgene run haplocheck@1.3.2 --files {params.input} --format vcf --output {params.out} --threads 1"

rule haplocheck_raw:
    input: rules.merge_vcfs.output
    output:
        'data/mity/results/haplocheck_res.txt',
        'data/mity/results/haplocheck_res.raw.txt',
        'data/mity/results/haplocheck_res.html'
    params:
        out = 'data/mity/results/haplocheck_res.txt'
    shell: "./src/haplocheck/haplocheck --out {params.out} {input} --raw"
