library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(purrr)
library(logger)

log_appender(appender_file(snakemake@log[[1]], append = TRUE))
log_threshold(TRACE)

coverage = snakemake@input[['coverage']]
vcf = snakemake@input[['vcf']]
samples = snakemake@params[['samples']]
mtdnacn.paths = snakemake@input[['mtdnacn']]
haplocheck.paths = snakemake@input[['haplocheck']]

## input tsv coverage file for annotate_coverage.py
log_info("Writing out coverage meterics file for annotate_coverage.py")
coverage_tsv = tibble(
  participant_id = samples,
  base_level_coverage_metrics = coverage,
  sample = samples
)

write_tsv(coverage_tsv, snakemake@output[['coverage_tsv']])

## participant vcf data file for combine_vcfs.py
log_info("Writing out participant vcf data file for combine_vcfs.py")
vcf_tsv = tibble(
  'entity:participant_id' = samples,
  s = samples,
  VCF = vcf,
)

write_tsv(vcf_tsv, snakemake@output[['vcf_tsv']])

## Sample Meta data
### entity:participant_id: Uploaded to Terra by user
### contamination: Output by Mutect2, gives the estimate of mitochondrial contamination
### freemix_percentage: Uploaded to Terra by user, can be calculated with VerifyBamID
### major_haplogroup: Output by Mutect2 which utilizes Haplogrep
### wgs_median_coverage: Uploaded to Terra by user, can be calculated with Picards CollectWgsMetrics
### mt_mean_coverage: Output by Mutect2, gives the mean mitochondrial coverage
log_info("Writing out Sample Meta data file for combine_vcfs.py")
mtdnacn <- map_df(mtdnacn.paths, read_tsv,
                  col_types = cols(
                    SampleID = col_character(),
                    model = col_character(),
                    mt_coverage = col_double(),
                    autosomal_coverage = col_double(),
                    mtcn_avg = col_double()
                  )) %>%
  pivot_wider(names_from = model, values_from = c('mt_coverage', 'autosomal_coverage', 'mtcn_avg'),) %>%
  select(SampleID, wgs_median_coverage = autosomal_coverage_median, mt_mean_coverage = mt_coverage_mean)

haplocheck <- map_df(haplocheck.paths, read_tsv, na = c("", "NA", "ND"),
                     col_types = cols(
                       Sample = col_character(),
                       `Contamination Status` = col_character(),
                       `Contamination Level` = col_double(),
                       Distance = col_double(),
                       `Sample Coverage` = col_double(),
                       `Overall Homoplasmies` = col_double(),
                       `Overall Heteroplasmies` = col_double(),
                       `Major Heteroplasmy Level` = col_double(),
                       `Minor Heteroplasmy Level` = col_double(),
                       `Major Haplogroup` = col_character(),
                       `Major Haplogroup Quality` = col_double(),
                       `Minor Haplogroup` = col_character(),
                       `Minor Haplogroup Quality` = col_double(),
                       `Major Homoplasmies Count` = col_double(),
                       `Minor Homoplasmies Count` = col_double(),
                       `Major Heteroplasmies Count` = col_double(),
                       `Minor Heteroplasmies Count` = col_double(),
                       Clusters = col_character()
                     )) %>%
  replace_na(list(`Contamination Level` = 0)) %>%
 select(Sample, contamination = `Contamination Level`, major_haplogroup = `Major Haplogroup`)

meta = left_join(haplocheck, mtdnacn, by = c("Sample" = "SampleID")) %>%
  mutate(freemix_percentage = 0) %>%
  select(s = Sample, 'entity:participant_id' = Sample, contamination, freemix_percentage, major_haplogroup, wgs_median_coverage, mt_mean_coverage)

write_tsv(meta, snakemake@output[['meta']])
