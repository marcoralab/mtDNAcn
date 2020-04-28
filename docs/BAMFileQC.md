# QC of AMP-AD BAM files

## ADNI
**Assembly of 809 whole mitochondrial genomes with clinical, imaging, and fluid biomarker phenotyping**

- ADNI was mapped to Hg19 which uses of the mitochondrial genome, represented as chrM, corresponding to NC_001807
- Since chrM and NC_012920 only differ by a few bases, we were able to extract only those reads that mapped to chrM (with SAMtools [46]), rather than all reads corresponding to the whole nuclear and mitochondrial genomes.
- Extracted reads were remapped to NC_012920 using BurrowsWheeler Aligner.
- performed local realignments around indels and base recalibration, which are not affected by ploidy, with Genome Analysis Toolkit to refine the new mappings.
- Used FreeBayes (-p 1 –F 0.6, and removed variants with quality less than 20) [47] to joint-call variants
- converted the resulting variant call format (VCF) file to fasta with vcf2fasta (vcflib, https://github.com/ vcflib/vcflib).
- annotated mitochondrial haplotypes with PhyMer [49]. Phy-Mer reports the five most likely mitochondrial haplotypes and a score, where 1 is a perfect score. For each of the samples, we selected the top hit. All samples had scores .0.99 except for one that had a score of 0.988.

## AMP-AD

[ROSMAP](https://www.synapse.org/#!Synapse:syn10901595), [Mayo](https://www.synapse.org/#!Synapse:syn10901601) and [MSBB](https://www.synapse.org/#!Synapse:syn10901600) WGS germline analysis
- Whole Genome data are processed on NYGC automated pipeline.
- Paired-end 150bp reads were aligned to the GRCh37 human reference using the Burrows-Wheeler Aligner (BWA-MEM v0.7.08)
- processed using the GATK best-practices workflow that includes marking of duplicate reads by the use of Picard tools v1.83, local realignment around indels, and base quality score recalibration (BQSR) via Genome Analysis Toolkit (GATK v3.4.0).

## DIAN
- Aligned to Hg19 / hg38
  - [Human Reference Discrepancies](https://gatk.broadinstitute.org/hc/en-us/articles/360035890711?id=23390#comparison)
  - GRCh38 and CRCh37 both use the rCRS as the reference for the mitochondrial genome

# QC of BAM Files

Runs the following:
1 - Picard CollectMultipleMetrics - output: ${sampleID}.metrics.insert_size_metrics
  - Collect multiple classes of metrics. This 'meta-metrics' tool runs one or more of the metrics collection modules at the same time to cut down on the time spent reading in data from input files. Available modules include CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution, MeanQualityByCycle, CollectBaseDistributionByCycle, CollectGcBiasMetrics, RnaSeqMetrics, CollectSequencingArtifactMetrics, and CollectQualityYieldMetrics. The tool produces outputs of '.pdf' and '.txt' files for each module, except for the CollectAlignmentSummaryMetrics module, which outputs only a '.txt' file. Output files are named by specifying a base name (without any file extensions).
  - ```java -jar $PICARD CollectMultipleMetrics I=temp/716.chrM.bam O=multiple_metrics```


2 - Samtools flagstat - output: ${sampleID}.flagstat
  – counts the number of alignments for each FLAG type
  - [flagstat](http://www.htslib.org/doc/samtools-flagstat.html)
  - ```samtools flagstat temp/716.chrM.bam > 716.flagstat.txt```

3 - Bamtools stats - output: ${sampleID}.stats
  - Description: prints general alignment statistics.
  - ```bamtools stats -in temp/716.chrM.bam```

4 - Picard CollectWgsMetrics - output: ${sampleID}.wgs
  - Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments. This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses. Both minimum base- and mapping-quality values as well as the maximum read depths (coverage cap) are user defined.
  - ```java -jar $PICARD CollectWgsMetrics I=temp/716.chrM.bam O=collect_wgs_metrics.txt R=raw/reference/Homo_sapiens_assembly19.fasta```

5 - Picard EstimateLibraryComplexity - output: ${sampleID}.complexity
  - Estimates the numbers of unique molecules in a sequencing library. This tool outputs quality metrics for a sequencing library preparation.Library complexity refers to the number of unique DNA fragments present in a given library. Reductions in complexity resulting from PCR amplification during library preparation will ultimately compromise downstream analyses via an elevation in the number of duplicate reads. PCR-associated duplication artifacts can result from: inadequate amounts of starting material (genomic DNA, cDNA, etc.), losses during cleanups, and size selection issues. Duplicate reads can also arise from optical duplicates resulting from sequencing-machine optical sensor artifacts.

6 - Sex Check - output: ${sampleID}.sexCheck
7 - Aneuploidy Check - output: ${sampleID}.aneuploidyCheck
