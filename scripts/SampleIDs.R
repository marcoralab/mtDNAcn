library(tidyverse)

setwd('/Users/sheaandrews/LOAD_minerva/dummy/shea/Projects/mDNACN')
files <- list.files('/Users/sheaandrews/LOAD_minerva/dummy/shea/Projects/mDNACN/raw') %>%
  tibble::enframe(name = NULL) %>%
  filter(!str_detect(value, '.crai')) %>%
  filter(!str_detect(value, '.bai')) %>%
  separate(value, c('sample', 'ext'), sep = '.final.') %>%
  arrange(ext) %>%
  filter(!is.na(ext))

files %>%
  select(sample) %>%
  write_tsv(., 'data/sampleIDs_test.txt' , col_names = F)

## Sample Info
## ROSMAP
rosmap.wgsqc <- read_csv('data/AMPAD_extra/rosmap/WGS_sample_QC_info.csv', guess_max = 10000)
rosmap.samp <- read_csv('data/AMPAD_extra/rosmap/ROSMAP_Clinical_2019-05_v3.csv') %>%
  left_join(rosmap.wgsqc, by = 'projid') %>%
  filter(!is.na(WGS_id)) %>%
  rename(study = Study) %>%
  select(SampleID = WGS_id, study)

## MSBB
msbb.samp <- read_tsv('data/AMPAD_extra/msbb.wgs.meta.tsv') %>%
  mutate(study = 'MSBB',
         Libid  = as.character(Libid )) %>%
  select(SampleID = Libid , study)

## Mayo
mayo.samp <- read_tsv('data/AMPAD/mayo/WGS_Metadata.txt') %>%
  mutate(WGS_Participant_ID = as.character(WGS_Participant_ID), study = 'MAYO')  %>%
  select(SampleID = WGS_Participant_ID, study)

## DIAN
dian.samp <- read_tsv('data/DIAN/completeMANIFEST') %>%
  filter(str_detect(file.name, 'cram$')) %>%
  mutate(wgs_id = str_remove(file.name, '.final.cram'),
         study = 'DIAN') %>%
  select(SampleID = wgs_id, study)


## All Samples
samples <- bind_rows(rosmap.samp, msbb.samp, mayo.samp, dian.samp)

SampleInfo <- list.files(path = 'raw', full.names = T, recursive = F) %>%
  tibble::enframe(name = NULL) %>%
  filter(str_detect(value, '.final.')) %>%
  mutate(SampleID = str_extract(value, "(?<=raw/).*(?=.final)"),
         file = ifelse(str_detect(value, "bam|cram"), 'wgs', NA),
         file = ifelse(str_detect(value, "bai|crai"), 'index', file)) %>%
  pivot_wider(values_from = value, names_from = file) %>%
  left_join(samples, by = 'SampleID') %>%
  mutate(ext = str_extract(wgs, "(?<=.final.).*"),
         hg = ifelse(study == 'DIAN', 38, 37))

write_csv(SampleInfo, 'data/SampleInfo.csv')  
