# snakejob -s MitochondriaPipeline.smk -j 10 --use-conda
# snakemake -s MitochondriaPipeline.smk --forceall --rulegraph | dot -Tpdf > MitochondriaPipeline.pdf

# Import Python Modules
import os
from itertools import product
import pandas as pd

SAMPLES = pd.read_csv('data/SampleInfo_test.csv', index_col='SampleID')
# ID = ["SM-CJGNJ", "SM-CJEIA", "71992"]
# ID = ["TWDF-037-0006-0106925-0010007709", "SM-CJGNJ"]

# Configfile
configfile: "config/gatk_config.yaml"
OUTDIR = config["outdir"]

## GATK Options
max_read_length = config["max_read_length"]
coverage_cap = config["coverage_cap"]
max_low_het_sites = config["max_low_het_sites"]
max_reads_per_alignment_start = config["max_reads_per_alignment_start"]
m2_extra_args = config["m2_extra_args"]
m2_extra_filtering_args = config["m2_extra_filtering_args"]
max_alt_allele_count = config["max_alt_allele_count"]
vaf_filter_threshold = config["vaf_filter_threshold"]
f_score_beta = config["f_score_beta"]
vaf_cutoff = config["vaf_cutoff"]

## Reference files
### hg38 whole genome reference
hg38_ref_fasta = config["reference"]["hg38"]["ref_fasta"]
hg38_ref_dict = config["reference"]["hg38"]["ref_dict"]
hg38_ref_fasta_index = config["reference"]["hg38"]["ref_fasta_index"]
### hg19 whole genome reference
hg19_ref_fasta = config["reference"]["hg19"]["ref_fasta"]
hg19_ref_dict = config["reference"]["hg19"]["ref_dict"]
hg19_ref_fasta_index = config["reference"]["hg19"]["ref_fasta_index"]
### hg38 mt refence files
mt_dict = config["reference"]["mt"]["mt_dict"]
mt_fasta = config["reference"]["mt"]["mt_fasta"]
mt_fasta_index = config["reference"]["mt"]["mt_fasta_index"]
mt_amb = config["reference"]["mt"]["mt_amb"]
mt_ann = config["reference"]["mt"]["mt_ann"]
mt_bwt = config["reference"]["mt"]["mt_bwt"]
mt_pac = config["reference"]["mt"]["mt_pac"]
mt_sa = config["reference"]["mt"]["mt_sa"]
non_control_region_interval_list = config["reference"]["mt"]["non_control_region_interval_list"]
blacklisted_sites = config["reference"]["mt"]["blacklisted_sites"]
blacklisted_sites_index = config["reference"]["mt"]["blacklisted_sites_index"]
### hg38 mt shifted reference files
mt_shifted_dict = config["reference"]["mt_shifted"]["mt_shifted_dict"]
mt_shifted_fasta = config["reference"]["mt_shifted"]["mt_shifted_fasta"]
mt_shifted_fasta_index = config["reference"]["mt_shifted"]["mt_shifted_fasta_index"]
mt_shifted_amb = config["reference"]["mt_shifted"]["mt_shifted_amb"]
mt_shifted_ann = config["reference"]["mt_shifted"]["mt_shifted_ann"]
mt_shifted_bwt = config["reference"]["mt_shifted"]["mt_shifted_bwt"]
mt_shifted_pac = config["reference"]["mt_shifted"]["mt_shifted_pac"]
mt_shifted_sa = config["reference"]["mt_shifted"]["mt_shifted_sa"]
shift_back_chain = config["reference"]["mt_shifted"]["shift_back_chain"]
control_region_shifted_reference_interval_list = config["reference"]["mt_shifted"]["control_region_shifted_reference_interval_list"]

# Varibles
RWD = os.getcwd()
BAM = ['bam', 'bai']

# GATK Variant Calling pipeline
## Adapated from https://github.com/gatk-workflows/gatk4-mitochondria-pipeline
## Laricchia, K. M. et al. (2021). BioRxiv. doi:10.1101/2021.07.23.453510

rule all:
    input:
        expand("{outdir}/{SAMPLE_IDS}/mosdepth/{SAMPLE_IDS}_mtdnacn.txt", outdir = OUTDIR, SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_per_base_coverage.tsv", outdir = OUTDIR, SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.final.split.vcf.gz", outdir = OUTDIR, SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("{outdir}/final/sample_vcf.vcf.bgz", outdir = OUTDIR)



def ref_fasta_input(wildcards):
    if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37:
        return hg19_ref_fasta
    elif SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38:
        return hg38_ref_fasta
    else:
        return print(SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] + ': Build Not Recognised')

def ref_index_input(wildcards):
    if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37:
        return hg19_ref_fasta_index
    elif SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38:
        return hg38_ref_fasta_index
    else:
        return print(SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] + ': Build Not Recognised')

def ref_dict_input(wildcards):
    if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37:
        return hg19_ref_dict
    elif SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38:
        return hg38_ref_dict
    else:
        return print(SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] + ': Build Not Recognised')

def mito_contig(wildcards):
    if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37:
        return 'MT'
    elif SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38:
        return 'chrM'
    else:
        return print(SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] + ': Build Not Recognised')

# rule CollectWgsMetrics_wgs:
#     input:
#         input_bam = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['wgs'],
#         input_bam_index = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['index'],
#         ref_fasta = ref_fasta_input,
#         ref_fasta_index = ref_index_input,
#     output:
#         meterics = '{outdir}/{SAMPLE_IDS}/{SAMPLE_IDS}_WGS_metrics.txt',
#         sens = '{outdir}/{SAMPLE_IDS}/{SAMPLE_IDS}_WGS_theoretical_sensitivity.txt'
#     params:
#         read_length = max_read_length,
#         coverage_cap = 100000,
#     conda: 'envs/gatk.yaml'
#     shell:
#         r"""
#         gatk CollectWgsMetrics \
#           -INPUT {input.input_bam} \
#           -VALIDATION_STRINGENCY SILENT \
#           -REFERENCE_SEQUENCE {input.ref_fasta} \
#           -OUTPUT {output.meterics} \
#           -USE_FAST_ALGORITHM true \
#           -READ_LENGTH {params.read_length} \
#           -COVERAGE_CAP {params.coverage_cap} \
#           -THEORETICAL_SENSITIVITY_OUTPUT WGS_theoretical_sensitivity.txt
#         """

rule mosdepth_median:
    input:
        input_bam = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['wgs'],
        input_bam_index = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['index'],
        ref_fasta = ref_fasta_input,
        ref_fasta_index = ref_index_input,
    output: "{outdir}/{SAMPLE_IDS}/mosdepth/{SAMPLE_IDS}.meadian.mosdepth.summary.txt"
    params:
        wd = '{outdir}/{SAMPLE_IDS}/mosdepth/{SAMPLE_IDS}.meadian'
    conda: 'envs/mosdepth.yaml'
    shell: "mosdepth -n --use-median -t 4 -f {input.ref_fasta} {params.wd} {input.input_bam}"

rule mosdepth_mean:
    input:
        input_bam = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['wgs'],
        input_bam_index = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['index'],
        ref_fasta = ref_fasta_input,
        ref_fasta_index = ref_index_input,
    output: "{outdir}/{SAMPLE_IDS}/mosdepth/{SAMPLE_IDS}.mean.mosdepth.summary.txt"
    params:
        wd = '{outdir}/{SAMPLE_IDS}/mosdepth/{SAMPLE_IDS}.mean'
    conda: 'envs/mosdepth.yaml'
    shell: "mosdepth -n -t 4 -f {input.ref_fasta} {params.wd} {input.input_bam}"

rule mtdnacn:
    input:
        median = rules.mosdepth_median.output,
        mean = rules.mosdepth_mean.output
    output:
        mtdnacn = "{outdir}/{SAMPLE_IDS}/mosdepth/{SAMPLE_IDS}_mtdnacn.txt"
    params:
        id = "{SAMPLE_IDS}"
    conda: 'envs/r.yaml'
    script: "scripts/mtdnacn.R"

# "Subsets a whole genome bam to just Mitochondria reads"
rule SubsetBamToChrM:
    input:
        input_bam = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['wgs'],
        input_bam_index = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['index'],
        ref_fasta = ref_fasta_input,
        ref_fasta_index = ref_index_input,
        ref_dict = ref_dict_input,
    params:
        contig_name = mito_contig,
    output:
        bam = '{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.bam',
        bai = '{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.bai',
    conda: 'envs/gatk.yaml'
    shell:
        r"""
        gatk PrintReads \
          -R {input.ref_fasta} \
          -L {params.contig_name} \
          --read-filter MateOnSameContigOrNoMappedMateReadFilter \
          --read-filter MateUnmappedAndUnmappedReadFilter \
          -I {input.input_bam} \
          --read-index {input.input_bam_index} \
          -O {output.bam}
        """

rule RevertSam:
    input:
        input_bam = rules.SubsetBamToChrM.output.bam,
        input_bai = rules.SubsetBamToChrM.output.bai
    output:
        unmapped_bam = temp('{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_unmapped.bam'),
    conda: 'envs/gatk.yaml'
    shell:
        r"""
        gatk RevertSam \
            -INPUT {input.input_bam} \
            -OUTPUT_BY_READGROUP false \
            -OUTPUT {output.unmapped_bam} \
            -VALIDATION_STRINGENCY LENIENT \
            -ATTRIBUTE_TO_CLEAR FT \
            -ATTRIBUTE_TO_CLEAR CO \
            -SORT_ORDER queryname \
            -RESTORE_ORIGINAL_QUALITIES false
        """

module AlignmentPipeline:
    snakefile: "AlignmentPipeline.smk"

use rule AlignAndMarkDuplicates from AlignmentPipeline as AlignToMt with:
    input:
        unmapped_bam = rules.RevertSam.output.unmapped_bam,
        mt_ref_fasta = mt_fasta
    output:
        fastq = temp('{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_unmapped.fastq'),
        aln_sam = temp('{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_aln.sam'),
        mba_bam = temp('{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_mba.bam'),
        md_bam = temp('{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_md.bam'),
        metrics_filename = '{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.duplicate_metrics',
        sorted_bam = '{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_sorted.bam',
        sorted_bai = '{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_sorted.bai'
    params:
        bwa_version = "0.7.17",
        bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 2 -Y",

use rule AlignAndMarkDuplicates from AlignmentPipeline as AlignToShiftedMt with:
    input:
        unmapped_bam = rules.RevertSam.output.unmapped_bam,
        mt_ref_fasta = mt_shifted_fasta
    output:
        fastq = temp('{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_shifted_unmapped.fastq'),
        aln_sam = temp('{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_shifted_aln.sam'),
        mba_bam = temp('{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_shifted_mba.bam'),
        md_bam = temp('{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_shifted_md.bam'),
        metrics_filename = '{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_shifted.duplicate_metrics',
        sorted_bam = '{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_shifted_sorted.bam',
        sorted_bai = '{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_shifted_sorted.bai'
    params:
        bwa_version = "0.7.17",
        bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 2 -Y",

rule CollectWgsMetrics:
    input:
        sorted_bam = rules.AlignToMt.output.sorted_bam,
        sorted_bam_index = rules.AlignToMt.output.sorted_bai,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
    output:
        meterics = '{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_metrics.txt',
        sens = '{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_theoretical_sensitivity.txt'
    params:
        read_length = max_read_length,
        coverage_cap = coverage_cap,
    conda: 'envs/gatk.yaml'
    shell:
        r"""
        gatk CollectWgsMetrics \
          -INPUT {input.sorted_bam} \
          -VALIDATION_STRINGENCY SILENT \
          -REFERENCE_SEQUENCE {input.mt_ref_fasta} \
          -OUTPUT {output.meterics} \
          -USE_FAST_ALGORITHM true \
          -READ_LENGTH {params.read_length} \
          -COVERAGE_CAP {params.coverage_cap} \
          -INCLUDE_BQ_HISTOGRAM true \
          -THEORETICAL_SENSITIVITY_OUTPUT {output.sens}
        """

rule CallMt:
    input:
        sorted_bam = rules.AlignToMt.output.sorted_bam,
        sorted_bai = rules.AlignToMt.output.sorted_bai,
        wgsmetrics = rules.CollectWgsMetrics.output.meterics,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
    output:
        vcf = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt.vcf.gz"),
        vcf_idx = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt.vcf.gz.tbi"),
        stats = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt.vcf.gz.stats"),
        bamOut = temp(touch("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt.bam"))
    params:
        max_reads_per_alignment_start = max_reads_per_alignment_start,
        m2_extra_args = m2_extra_args,
        region = "chrM:576-16024"
    conda: 'envs/gatk.yaml'
    shell:
        r"""
        gatk Mutect2 \
            -R {input.mt_ref_fasta} \
            -I {input.sorted_bam} \
            --read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --read-filter MateUnmappedAndUnmappedReadFilter \
            -O {output.vcf} \
            --bam-output {output.bamOut} \
            {params.m2_extra_args} -L {params.region} \
            --annotation StrandBiasBySample \
            --mitochondria-mode \
            --max-reads-per-alignment-start {params.max_reads_per_alignment_start} \
            --max-mnp-distance 0
        """

use rule CallMt as CallShiftedMt with:
    input:
        sorted_bam = rules.AlignToShiftedMt.output.sorted_bam,
        sorted_bai = rules.AlignToShiftedMt.output.sorted_bai,
        wgsmetrics = rules.CollectWgsMetrics.output.meterics,
        mt_ref_fasta = mt_shifted_fasta,
        mt_ref_fasta_index = mt_shifted_fasta_index,
        mt_dict = mt_shifted_dict,
    output:
        vcf = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt_shifted.vcf.gz"),
        vcf_idx = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt_shifted.vcf.gz.tbi"),
        stats = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt_shifted.vcf.gz.stats"),
        bamOut = temp(touch("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt_shifted.bam")),
    params:
        max_reads_per_alignment_start = max_reads_per_alignment_start,
        m2_extra_args = m2_extra_args,
        region = "chrM:8025-9144"


rule LiftoverAndCombineVcfs:
    input:
        shifted_vcf = rules.CallShiftedMt.output.vcf,
        vcf = rules.CallMt.output.vcf,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
        shift_back_chain = shift_back_chain
    output:
        shifted_back_vcf = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}.shifted_back.vcf"),
        rejected_vcf = "{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}.rejected.vcf",
        merged_vcf = "{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}.merged.vcf",
        merged_vcf_index = "{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}.merged.vcf.idx"
    conda: 'envs/gatk.yaml'
    shell:
        r"""
        gatk LiftoverVcf \
          -I {input.shifted_vcf} \
          -O {output.shifted_back_vcf} \
          -R {input.mt_ref_fasta} \
          -CHAIN {input.shift_back_chain} \
          -REJECT {output.rejected_vcf}

        gatk MergeVcfs \
          -I {output.shifted_back_vcf} \
          -I {input.vcf} \
          -O {output.merged_vcf}
        """

rule MergeStats:
    input:
        shifted_stats = rules.CallMt.output.stats,
        non_shifted_stats = rules.CallShiftedMt.output.stats
    output:
        "{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_raw_combined.stats"
    conda: 'envs/gatk.yaml'
    shell:
        r"""
        gatk MergeMutectStats --stats {input.shifted_stats} --stats {input.non_shifted_stats} -O {output}
        """

rule InitialFilter:
    input:
        raw_vcf = rules.LiftoverAndCombineVcfs.output.merged_vcf,
        raw_vcf_index = rules.LiftoverAndCombineVcfs.output.merged_vcf_index,
        raw_vcf_stats = rules.MergeStats.output,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
        blacklisted_sites = blacklisted_sites,
        blacklisted_sites_index = blacklisted_sites_index
    output:
        FilterMutectCalls_vcf = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterMutectCalls.vcf"),
        filtered_vcf = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_InitialFilter.vcf.gz"),
        filtered_vcf_idx = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_InitialFilter.vcf.gz.tbi"),
        bamout = temp(touch("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_InitialFilter.bam")) # We need to create these files regardless, even if they stay empty
    params:
        m2_extra_filtering_args = m2_extra_filtering_args,
        max_alt_allele_count = max_alt_allele_count,
        vaf_filter_threshold = vaf_filter_threshold,
        f_score_beta = f_score_beta,
    conda: 'envs/gatk.yaml'
    shell:
        r"""
        gatk  FilterMutectCalls \
            -V {input.raw_vcf} \
            -R {input.mt_ref_fasta} \
            -O {output.FilterMutectCalls_vcf} \
            --stats {input.raw_vcf_stats} \
            {params.m2_extra_filtering_args} \
            --max-alt-allele-count {params.max_alt_allele_count} \
            --mitochondria-mode \
            --min-allele-fraction {params.vaf_filter_threshold} \
            --f-score-beta {params.f_score_beta}

        gatk VariantFiltration \
            -V {output.FilterMutectCalls_vcf} \
            -O {output.filtered_vcf} \
            --apply-allele-specific-filters \
            --mask {input.blacklisted_sites} \
            --mask-name "blacklisted_site"
        """

rule SplitMultiAllelicsAndRemoveNonPassSites:
    input:
        filtered_vcf = rules.InitialFilter.output.filtered_vcf,
        filtered_vcf_idx = rules.InitialFilter.output.filtered_vcf_idx,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
    output:
        split_vcf = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_split.vcf"),
        vcf_for_haplochecker = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_splitAndPassOnly.vcf"),
    conda: 'envs/gatk.yaml'
    shell:
        r"""
        gatk LeftAlignAndTrimVariants \
          -R {input.mt_ref_fasta} \
          -V {input.filtered_vcf} \
          -O {output.split_vcf} \
          --split-multi-allelics \
          --dont-trim-alleles \
          --keep-original-ac

          gatk SelectVariants \
            -V {output.split_vcf} \
            -O {output.vcf_for_haplochecker} \
            --exclude-filtered
        """

rule GetContamination:
    input:
        filtered_vcf = rules.SplitMultiAllelicsAndRemoveNonPassSites.output.vcf_for_haplochecker,
    output:
        contamination_tmp = temp("{outdir}/{SAMPLE_IDS}/haplocheck/{SAMPLE_IDS}_contamination.txt"),
        contamination_raw_tmp = temp("{outdir}/{SAMPLE_IDS}/haplocheck/{SAMPLE_IDS}_contamination.raw.txt"),
        contamination = "{outdir}/{SAMPLE_IDS}/haplocheck/{SAMPLE_IDS}_chrM_contamination.txt",
        contamination_raw = "{outdir}/{SAMPLE_IDS}/haplocheck/{SAMPLE_IDS}_chrM_contamination.raw.txt",
    shell:
        r"""
        java -jar src/haplocheck/haplocheck.jar --raw --out {output.contamination_tmp} {input}
        sed 's/\"//g' {output.contamination_raw_tmp} > {output.contamination_raw}
        sed 's/\"//g' {output.contamination_tmp} > {output.contamination}
        """

rule FilterContamination:
    input:
        filtered_vcf = rules.InitialFilter.output.filtered_vcf,
        filtered_vcf_idx = rules.InitialFilter.output.filtered_vcf_idx,
        raw_vcf_stats = rules.MergeStats.output,
        contamination_raw = rules.GetContamination.output.contamination_raw,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
        blacklisted_sites = blacklisted_sites,
        blacklisted_sites_index = blacklisted_sites_index,
    output:
        FilterContaminationMutectCalls_vcf = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterContaminationMutectCalls.vcf"),
        filtered_vcf = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterContamination.vcf.gz"),
        filtered_vcf_idx = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterContamination.vcf.gz.tbi"),
        bamout = temp(touch("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterContamination.bam"))
    params:
        f_score_beta = f_score_beta,
        m2_extra_filtering_args = m2_extra_filtering_args,
        max_alt_allele_count = max_alt_allele_count,
        vaf_filter_threshold = vaf_filter_threshold,
    conda: 'envs/gatk.yaml'
    shell:
        r"""
        hc_contamination=$(awk -F "\t" 'FNR == 2' {input.contamination_raw} | awk '/YES/{{print $3;next}}{{print "0"}}')

        gatk  FilterMutectCalls \
            -V {input.filtered_vcf} \
            -R {input.mt_ref_fasta} \
            -O {output.FilterContaminationMutectCalls_vcf} \
            --stats {input.raw_vcf_stats} \
            {params.m2_extra_filtering_args} \
            --max-alt-allele-count {params.max_alt_allele_count} \
            --mitochondria-mode \
            --f-score-beta {params.f_score_beta}
            --contamination-estimate $hc_contamination

        gatk VariantFiltration \
            -V {output.FilterContaminationMutectCalls_vcf} \
            -O {output.filtered_vcf} \
            --apply-allele-specific-filters \
            --mask {input.blacklisted_sites} \
            --mask-name "blacklisted_site"
        """

rule FilterNuMTs:
    input:
        filtered_vcf = rules.FilterContamination.output.filtered_vcf,
        filtered_vcf_idx = rules.FilterContamination.output.filtered_vcf_idx,
        autosomal_coverage = rules.mtdnacn.output,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
    output:
        numt_filtered_vcf = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterNuMTs.vcf.gz"),
        numt_filtered_vcf_idx = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterNuMTs.vcf.gz.tbi"),
    conda: 'envs/gatk.yaml'
    shell:
        r"""
        auto_cov=$(awk '{{ if ($2 == "median") {{ print $4}} }}' {input.autosomal_coverage})
        echo $auto_cov
        gatk NuMTFilterTool \
          -R {input.mt_ref_fasta} \
          -V {input.filtered_vcf} \
          -O {output.numt_filtered_vcf} \
          --autosomal-coverage $auto_cov
        """

rule FilterLowHetSites:
    input:
        filtered_vcf = rules.FilterNuMTs.output.numt_filtered_vcf,
        filtered_vcf_idx = rules.FilterNuMTs.output.numt_filtered_vcf_idx,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
    output:
        final_filtered_vcf = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.final.vcf.gz"),
        final_filtered_vcf_idx = temp("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.final.vcf.gz.tbi"),
    params:
        max_sites = max_low_het_sites,
    conda: 'envs/gatk.yaml'
    shell:
        r"""
        gatk MTLowHeteroplasmyFilterTool \
          -R {input.mt_ref_fasta} \
          -V {input.filtered_vcf} \
          -O {output.final_filtered_vcf} \
          --max-allowed-low-hets {params.max_sites}
        """

rule SplitMultiAllelicSites:
    input:
        filtered_vcf = rules.FilterLowHetSites.output.final_filtered_vcf,
        filtered_vcf_idx = rules.FilterLowHetSites.output.final_filtered_vcf_idx,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
    output:
        split_vcf = "{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.final.split.vcf.gz",
        split_vcf_idx = "{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.final.split.vcf.gz.tbi",
    conda: 'envs/gatk.yaml'
    shell:
        r"""
        gatk LeftAlignAndTrimVariants \
          -R {input.mt_ref_fasta} \
          -V {input.filtered_vcf} \
          -O {output.split_vcf} \
          --split-multi-allelics \
          --dont-trim-alleles \
          --keep-original-ac
        """

rule CoverageAtEveryBase:
    input:
        bam_regular_ref = rules.AlignToMt.output.sorted_bam,
        bam_regular_ref_index = rules.AlignToMt.output.sorted_bai,
        bam_shifted_ref = rules.AlignToShiftedMt.output.sorted_bam,
        bam_shifted_ref_index = rules.AlignToShiftedMt.output.sorted_bai,
        shift_back_chain = shift_back_chain,
        control_region_shifted_reference_interval_list = control_region_shifted_reference_interval_list,
        non_control_region_interval_list = non_control_region_interval_list,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
        shifted_mt_ref_fasta = mt_shifted_fasta,
        shifted_mt_ref_fasta_index = mt_shifted_fasta_index,
        shifted_mt_dict = mt_shifted_dict,
    output:
        table = "{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_per_base_coverage.tsv",
        pbc_non_control_region_tsv = "{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_non_control_region.tsv",
        pbc_non_control_region_meterics = "{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_non_control_region.metrics",
        pbc_control_region_shifted_tsv = "{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_control_region_shifted.tsv",
        pbc_control_region_shifted_meterics = "{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_control_region_shifted.metrics"
    conda: 'envs/gatk.yaml'
    shell:
        r"""
        gatk CollectHsMetrics \
          -I {input.bam_regular_ref} \
          -R {input.mt_ref_fasta} \
          -PER_BASE_COVERAGE {output.pbc_non_control_region_tsv} \
          -O {output.pbc_non_control_region_meterics}\
          -TI {input.non_control_region_interval_list} \
          -BI {input.non_control_region_interval_list} \
          -covMax 20000 \
          -SAMPLE_SIZE 1

        gatk CollectHsMetrics \
          -I {input.bam_shifted_ref} \
          -R {input.shifted_mt_ref_fasta} \
          -PER_BASE_COVERAGE {output.pbc_control_region_shifted_tsv} \
          -O {output.pbc_control_region_shifted_meterics} \
          -TI {input.control_region_shifted_reference_interval_list} \
          -BI {input.control_region_shifted_reference_interval_list} \
          -covMax 20000 \
          -SAMPLE_SIZE 1

        Rscript workflow/scripts/CoverageAtEveryBase.R {output.pbc_non_control_region_tsv} {output.pbc_control_region_shifted_tsv} {output.table}

        """

rule hail_inputs:
    input:
        coverage = expand("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_per_base_coverage.tsv", outdir = OUTDIR, SAMPLE_IDS = SAMPLES.index.tolist()),
        vcf = expand("{outdir}/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.final.split.vcf.gz", outdir = OUTDIR, SAMPLE_IDS = SAMPLES.index.tolist()),
        mtdnacn = expand("{outdir}/{SAMPLE_IDS}/mosdepth/{SAMPLE_IDS}_mtdnacn.txt", outdir = OUTDIR, SAMPLE_IDS = SAMPLES.index.tolist()),
        haplocheck = expand("{outdir}/{SAMPLE_IDS}/haplocheck/{SAMPLE_IDS}_chrM_contamination.raw.txt", outdir = OUTDIR, SAMPLE_IDS = SAMPLES.index.tolist())
    output:
        coverage_tsv = "{outdir}/sample_coverage_paths.txt",
        vcf_tsv = "{outdir}/sample_vcf_paths.txt",
        meta = "{outdir}/sample_metadata.txt",
    params:
        samples = expand("{SAMPLE_IDS}", SAMPLE_IDS = SAMPLES.index.tolist())
    conda: "envs/r.yaml"
    script: "scripts/hail_inputs.R"


rule annotate_coverage:
    input:
        input_tsv = rules.hail_inputs.output.coverage_tsv,
        script = 'workflow/scripts/annotate_coverage.py'
    output:
        ht = directory("{outdir}/annotate_coverage/all_chrM_annotate_coverage.ht"),
        mt = directory("{outdir}/annotate_coverage/all_chrM_annotate_coverage.mt"),
        all = "{outdir}/annotate_coverage/all_chrM_annotate_coverage.tsv",
        sample = "{outdir}/annotate_coverage/all_chrM_annotate_coverage_sample_level.txt"
    params:
        temp_dir = "{outdir}/temp"
    conda: 'envs/hail.yaml'
    shell: "{input.script} -i {input.input_tsv} -o {output.ht} -t {params.temp_dir} --overwrite"

rule combine_vcfs:
    input:
        script = 'workflow/scripts/combine_vcfs.py',
        vcf_paths = rules.hail_inputs.output.vcf_tsv,
        coverage = rules.annotate_coverage.output.mt,
        ARTIFACT_PRONE_SITES = "raw/gnomad_resources/artifact_prone_sites.bed"
    output:
        vcf = "{outdir}/combine_vcfs/all_chrM.vcf.bgz",
        mt = directory("{outdir}/combine_vcfs/all_chrM.mt"),
        temp = temp(directory("{outdir}/combine_vcfs/temp"))
    params:
        output = "all_chrM",
        VCF_COL_NAME = "VCF",
        temp_dir = "{outdir}/combine_vcfs/temp",
        output_bucket = "{outdir}/combine_vcfs"
    conda: 'envs/hail.yaml'
    shell:
        r"""
        {input.script} --overwrite -p {input.vcf_paths} -c {input.coverage} \
        -v {params.VCF_COL_NAME} -a {input.ARTIFACT_PRONE_SITES} \
        -t {params.temp_dir} -f {params.output} -o {params.output_bucket}
        """

rule add_annotations:
    input:
        script = "workflow/scripts/add_annotations.py",
        vcf = rules.combine_vcfs.output.vcf,
        mt = rules.combine_vcfs.output.mt,
        meta = rules.hail_inputs.output.meta,
    output:
        mt = directory("{outdir}/final/annotated_combined.mt"),
        ht = directory("{outdir}/final/combined_sites_only.ht"),
        txt = "{outdir}/final/combined_sites_only.txt",
        sites_vcf = "{outdir}/final/combined_sites_only.vcf.bgz",
        sample = "{outdir}/final/sample_annotations.txt",
        sample_vcf = "{outdir}/final/sample_vcf.vcf.bgz",
    params:
        temp_dir = "{outdir}/temp",
        OUTPUT_DIR = "{outdir}/final",
        min_het_threshold = 0.10,
        min_hom_threshold = 0.95,
        vaf_filter_threshold = vaf_filter_threshold
    conda: 'envs/hail.yaml'
    shell:
        r"""
        {input.script} -m {input.mt} -d {params.OUTPUT_DIR} -a {input.meta} -v {params.temp_dir}\
        --min_het_threshold {params.min_het_threshold} --min_hom_threshold {params.min_hom_threshold} --vaf_filter_threshold {params.vaf_filter_threshold} \
        --run_vep --overwrite
        """
