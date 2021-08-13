library(tidyverse)
library(glue)
# ln -s /sc/arion/projects/AMPADWGS/RawData/*/*/*/*final.bam .
# ln -s /sc/arion/projects/AMPADWGS/RawData/*/*/*/*final.bai .
# ln -s /sc/arion/projects/AMPADWGS/RawDataMayo/*/*/*/*final.bam .
# ln -s /sc/arion/projects/AMPADWGS/RawDataMayo/*/*/*/*final.bai .
# ln -s /sc/arion/projects/AMPADWGS/RawDataSinai/*/*/*/*final.bam .
# ln -s /sc/arion/projects/AMPADWGS/RawDataSinai/*/*/*/*final.bai .
# ln -s /sc/arion/projects/LOAD/Data/DIAN/DIAN_WGS/xfer.genome.wustl.edu/gxfer1/40905687920482/*/*final.cram .
# ln -s /sc/arion/projects/LOAD/Data/DIAN/DIAN_WGS/xfer.genome.wustl.edu/gxfer1/40905687920482/*/*final.crai .

setwd('/sc/arion/projects/LOAD/shea/Projects/mtDNAcn')
files <- list.files('/sc/arion/projects/LOAD/shea/Projects/mtDNAcn/raw/seq') %>%
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
SampleInfo <- list.files(path = 'raw/seq', full.names = T, recursive = F) %>%
  tibble::enframe(name = NULL) %>%
  filter(str_detect(value, '.final.')) %>%
  mutate(SampleID = str_extract(value, "(?<=raw/seq/).*(?=.final)"),
         file = ifelse(str_detect(value, "bam|cram"), 'wgs', NA),
         file = ifelse(str_detect(value, "bai|crai"), 'index', file)) %>%
  pivot_wider(values_from = value, names_from = file)

## ROSMAP
rosmap.wgsqc <- read_csv('data/AMPAD_extra/rosmap/WGS_sample_QC_info.csv', guess_max = 10000)
rosmap.samp <- read_csv('data/AMPAD_extra/rosmap/ROSMAP_Clinical_2019-05_v3.csv') %>%
  left_join(rosmap.wgsqc, by = 'projid') %>%
  filter(!is.na(WGS_id)) %>%
  filter(!is.na(Study)) %>%
  rename(study = Study) %>%
  mutate(hg = 37) %>%
  select(SampleID = WGS_id, study, hg) %>%
  left_join(SampleInfo, by = 'SampleID')

## MSBB
msbb.samp <- read_tsv('data/AMPAD_extra/msbb.wgs.meta.tsv') %>%
  mutate(study = 'MSBB',
         Libid  = as.character(Libid ), 
         hg = 37) %>%
  select(SampleID = Libid , study, hg) %>%
  left_join(SampleInfo, by = 'SampleID') 

## Mayo
mayo.samp <- read_tsv('data/AMPAD/mayo/WGS_Metadata.txt') %>%
  mutate(WGS_Participant_ID = as.character(WGS_Participant_ID), 
         study = 'MAYO', 
         hg = 37)  %>%
  select(SampleID = WGS_Participant_ID, study, hg) %>%
  left_join(SampleInfo, by = 'SampleID')

## DIAN
dian.samp <- read_tsv('data/DIAN/completeMANIFEST') %>%
  filter(str_detect(file.name, 'cram$')) %>%
  mutate(wgs_id = str_remove(file.name, '.final.cram'),
         study = 'DIAN', 
         hg = 38) %>%
 select(SampleID = subject.name, wgs_id, study, hg) %>%
 left_join(SampleInfo, by = c('wgs_id' = 'SampleID')) %>%
  # select(SampleID =  wgs_id, study) %>%
  # left_join(SampleInfo, by = 'SampleID') %>%
  select(-wgs_id)

## All Samples
samples <- bind_rows(rosmap.samp, msbb.samp, mayo.samp, dian.samp) %>%
  mutate(
    ext = str_extract(wgs, "(?<=.final.).*"),
    # SampleID = glue("{study}-{SampleID}")
    )

write_csv(samples, 'data/SampleInfo.csv') 

## Testing dataset 
samples %>% group_by(study) %>% slice(1) %>% write_csv("sandbox/SampleInfo_test.csv")

samples %>% 
  filter(study %in% c("ROS", "MAP")) %>%
  write_csv(., 'data/SampleInfo_rosmap.csv') 

samples %>% 
  filter(study == "MAYO") %>%
  write_csv(., 'data/SampleInfo_mayo.csv') 

samples %>% 
  filter(study == "MSBB") %>%
  write_csv(., 'data/SampleInfo_msbb.csv') 

samples %>% 
  filter(study == "DIAN") %>%
  write_csv(., 'data/SampleInfo_dian.csv') 
