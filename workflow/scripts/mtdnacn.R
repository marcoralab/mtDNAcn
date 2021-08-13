#!/usr/bin/Rscript
# ============================================================================ #
#                     Estimate mtDNAcn from Mosdepth output
# Mosdepth is used to calcuate genome-wide sequencing coverage. mtDNAcn is
# estimated as the ratio of the average mitochodnrial DNA coverage by the 
# average autosomal DNA 
# ============================================================================ #

# Setup ------------------------------------------------------------------------
message("\n\nSetup for estimating mtDNAcn\n\n")
library(tidyverse)
library(glue)

mean.path = snakemake@input[["mean"]]
median.path = snakemake@input[["median"]]
sampleID = snakemake@params[["id"]]
outfile = snakemake@output[["mtdnacn"]]
 
# Sandbox ----------------------------------------------------------------------
# sampleID = "TWDF-035-0038-0300113-1161310248"
# mean.path = glue("sandbox/output/{sampleID}/mosdepth/{sampleID}.mean.mosdepth.summary.txt")
# median.path = glue("sandbox/output/{sampleID}/mosdepth/{sampleID}.meadian.mosdepth.summary.txt")
# outfile = snakemake@output[["mtdnacn"]]


# Load Data --------------------------------------------------------------------
message(glue::glue(
 "\n\nReading in file:  {mean.path} & {median.path}  \n\n"
  ))

men.raw <- read_tsv(mean.path)
med.raw <- read_tsv(median.path)

# Estimate mtDNAcn -------------------------------------------------------------
message("\n\nEstimating mtDNAcn\n\n")
dat.men <- men.raw %>%
  mutate(chrom = str_replace(chrom, 'chrM', 'MT'),
         chrom = str_replace(chrom, 'chr', ''),
         autosome = chrom %in% 1:22,
         SampleID = sampleID, 
         model = "mean") %>%
  filter(chrom %in% c(1:22, "MT")) %>%
  group_by(autosome) %>%
  summarise(SampleID = first(SampleID),
            model = first(model),
            coverage = mean(mean)) %>%
  mutate(autosome = case_when(autosome == FALSE ~ "mt_coverage",
                              autosome == TRUE ~ "autosomal_coverage")) %>%
  pivot_wider(names_from = autosome, values_from = coverage) %>%
  mutate(mtcn_avg = (mt_coverage / autosomal_coverage)*2)

dat.med <- med.raw %>%
  mutate(chrom = str_replace(chrom, 'chrM', 'MT'),
         chrom = str_replace(chrom, 'chr', ''),
         autosome = chrom %in% 1:22,
         SampleID = sampleID, 
         model = "median") %>%
  filter(chrom %in% c(1:22, "MT")) %>%
  group_by(autosome) %>%
  summarise(SampleID = first(SampleID),
            model = first(model),
            coverage = median(mean)) %>%
  mutate(autosome = case_when(autosome == FALSE ~ "mt_coverage",
                              autosome == TRUE ~ "autosomal_coverage")) %>%
  pivot_wider(names_from = autosome, values_from = coverage) %>%
  mutate(mtcn_avg = (mt_coverage / autosomal_coverage)*2)

dat <- bind_rows(dat.men, dat.med)
print(dat)

# Export -----------------------------------------------------------------------
message(glue::glue(
  "\n\nExporting file to: ", {outfile},  "\n\n"
))
write_tsv(dat, outfile)


























