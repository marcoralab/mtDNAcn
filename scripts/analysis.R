library(readxl)
library(forcats)
library(tidyverse)


setwd('/Users/sheaandrews/LOAD_minerva/dummy/shea/Projects/mDNACN')

mtdna <- list.files(path = "data/mitocalc", pattern = ".txt", full.names = TRUE) %>%
  map(., function(x){
    file <- x
    d <- read_delim(file, delim = ':  ', col_names = F) %>%
      mutate(X2 = str_trim(X2, side = "both"))
    tibble(sample = str_remove(str_remove(file, 'data/mitocalc/'), 'MTData.txt'),
           mtcn_avg = filter(d, X1 == 'mt_copy_number_avg') %>% pull(X2),
           mt_coverage = filter(d, X1 == 'mt_coverage') %>% pull(X2),
           autosomal_coverage = filter(d, X1 == 'autosomal_coverage') %>% pull(X2),
           actual_basepairs_covered_by_reads = filter(d, X1 == 'actual basepairs covered by reads') %>% pull(X2),
           hrom_used_for_autosomal_coverage = filter(d, X1 == 'chrom_used_for_autosomal_coverage') %>% pull(X2))

  }) %>%
  bind_rows()

## ROS/MAP
rosmap <- read_tsv('data/AMPAD/rosmap/WGS_Metadata.txt')

## Mount Sinai Brain Bank
msbb <- read_tsv('data/AMPAD/msbb/WGS_Metadata.txt') %>%
  unite(apoe, c('Apo1', 'Apo2'), sep = "")
msbb$dx = 'OTHER'
msbb$dx[msbb$CDR <= 0.5 & msbb$bbscore <= 3 & msbb$NP.1 <= 1] = 'CONTROL'
msbb$dx[msbb$CDR >= 1 & msbb$bbscore >= 4 & msbb$NP.1 >= 2] = 'AD'

## Mayo Clinic
mayo <- read_tsv('data/AMPAD/mayo/WGS_Metadata.txt')



rosmpaid <- read_csv('data/AMPAD/rosmap/fromSynapse/ROSMAP_IDkey.csv')
dat <- read_xls("~/Dropbox/Research/PhD/Analysis/10 CAPA/RUSH Cross Sectional 04-2015.xls")
rosmap_dx <- read_xls("~/Dropbox/Research/PhD/Analysis/10 CAPA/RUSH Longitudinal 04-2015.xls") %>%
  mutate(projid = str_remove(projid, "^0+")) %>%
  mutate(projid = as.numeric(projid)) %>%
  group_by(projid) %>%
  arrange(fu_year) %>%
  slice(which.max(fu_year)) %>%
  ungroup()

test <- rosmap %>%
  select(wgs_id, projid, msex, race, spanish, apoe_genotype, age_at_visit_max, cogdx) %>%
  arrange(projid) %>%
  left_join(rosmap_dx, by = 'projid') %>%
  left_join(select(mtdna, sample, mtcn_avg), by = c('wgs_id' = 'sample')) %>%
  mutate(mtcn_avg = as.numeric(mtcn_avg)) %>%
  mutate(study = ifelse(is.na(study), 'ROS', study)) %>%
  mutate(study = as.factor(study))

lm(age_at_visit ~ mtcn_avg + msex, data = test) %>%
  summary(.)
