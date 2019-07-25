## Snakefile for calculating mtDNA copy number

configfile: "config.yaml"
BAM_FILES = config['BamFile']

rule all:
    input:
        expand("data/{BAM_FILES}Data.txt", bam_id = BAM_FILES)

rule fastMitoCalc:
    input:
        bam = 'raw/{BAM_FILES}.final.bam',
        bai = 'raw/{BAM_FILES}.final.bai'
    output: "data/{BAM_FILES}Data.txt"
    shell:
        'perl src/fastMitoCalc/fastMitoCalc.pl -f {input.bam} -w data -p src/fastMitoCalc/BaseCoverage'
