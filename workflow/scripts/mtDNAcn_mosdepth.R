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
infile = snakemake@input[["mosdepth"]]
sampleID = snakemake@params[["id"]]
outfile = snakemake@output[["mtdnacn"]]

# Sandbox ----------------------------------------------------------------------
# setwd("/sc/arion/projects/LOAD/shea/Projects/mtDNAcn/")
# infile = "data/mosdepth/SM-CJEIA.mosdepth.summary.txt"
# sampleID = "SM-CJEIA"
# outfile =

# Load Data --------------------------------------------------------------------
message(glue::glue(
 "\n\nReading in file: ", {infile},  "\n\n"
  ))
dat.raw <- read_tsv(infile)

# Estimate mtDNAcn -------------------------------------------------------------
message("\n\nEstimating mtDNAcn\n\n")
dat <- dat.raw %>%
  mutate(chrom = str_replace(chrom, 'chrM', 'MT'),
         chrom = str_replace(chrom, 'chr', ''),
         autosome = chrom %in% 1:22,
         SampleID = sampleID) %>%
  filter(chrom %in% c(1:22, "MT")) %>%
  group_by(autosome) %>%
  summarise(SampleID = first(SampleID),
            coverage = mean(mean)) %>%
  mutate(autosome = case_when(autosome == FALSE ~ "mt_coverage",
                              autosome == TRUE ~ "autosomal_coverage")) %>%
  pivot_wider(names_from = autosome, values_from = coverage) %>%
  mutate(mtcn_avg = (mt_coverage / autosomal_coverage)*2)

print(dat)

# Export -----------------------------------------------------------------------
message(glue::glue(
  "\n\nExporting file to: ", {outfile},  "\n\n"
))
write_tsv(dat, outfile)


























