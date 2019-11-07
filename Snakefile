## Snakefile for calculating mtDNA copy number

configfile: "config.yaml"
# BAM_FILES = config['BamFile']
BAM_FILES = [line.rstrip('\n') for line in open('data/sampleIDs.txt')]

shell.prefix('module load samtools; ')

rule all:
    input:
        expand("data/mitocalc/{BAM_FILES}MTData.txt", BAM_FILES = BAM_FILES)

rule fastMitoCalc:
    input:
        bam = 'raw/{BAM_FILES}.final.bam',
        bai = 'raw/{BAM_FILES}.final.bai'
    output: "data/mitocalc/{BAM_FILES}MTData.txt"
    shell:
        'perl src/fastMitoCalc/fastMitoCalc.pl -f {input.bam} -w data/mitocalc -p src/fastMitoCalc/BaseCoverage'
