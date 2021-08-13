#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
args = commandArgs(trailingOnly=TRUE)

infile = args[1]

dat.raw <- suppressMessages(read_tsv(infile))
dat.raw %>%
  mutate(chrom = str_replace(chrom, 'chrM', 'MT'),
         chrom = str_replace(chrom, 'chr', '')) %>%
  filter(chrom %in% 1:22) %>%
  summarise(mean = mean(mean)) %>% pull(mean) %>% round(2) %>% cat()
