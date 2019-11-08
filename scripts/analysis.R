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
rosmap.raw <- read_tsv('data/AMPAD/rosmap/WGS_Metadata.txt') %>%
  left_join(mtdna, by = c('wgs_id' = 'sample'))

rosmap <- rosmap.raw %>%
  filter(race == 1) %>%
  mutate(aod = str_replace(age_death, '\\+', ''),
        aod = round(as.numeric(aod), 0),
        dx = recode(cogdx, '1' = 'CTRL', '4' = 'AD', '5' = 'AD'),
        dx = fct_relevel(dx, 'CTRL', 'AD'),
        apoe4 = recode(apoe_genotype, '22' = 'e4-', '23' = 'e4-', '33' = 'e4-', '24' = 'e4+', '34' = 'e4+', '44' = 'e4+'),
        aod_cat = cut(aod, c(50, 60, 70, 80, 90, Inf), c('50-69', '60-69', '70-79', '80-89', '90+'), right = FALSE),
        msex = as.factor(msex), dx = as.factor(dx), mtcn_avg = as.numeric(mtcn_avg), cogdx = as.factor(cogdx), apoe_genotype = as.factor(apoe_genotype), apoe4 = as.factor(apoe4)) %>%
  select(id = wgs_id, study, sex = msex, apoe = apoe_genotype, apoe4, aod, aod_cat, mtcn_avg, dx.raw = cogdx, dx)

## Mount Sinai Brain Bank
msbb.raw <- read_tsv('data/AMPAD/msbb/WGS_Metadata.txt') %>%
  unite(apoe, c('Apo1', 'Apo2'), sep = "") %>%
  mutate(WGS = as.character(WGS)) %>%
  left_join(mtdna, by = c('WGS' = 'sample'))

msbb.raw$dx.raw = 'OTHER'
msbb.raw$dx.raw[msbb.raw$CDR <= 0.5 & msbb.raw$bbscore <= 3 & msbb.raw$NP.1 <= 1] = 'CONTROL'
msbb.raw$dx.raw[msbb.raw$CDR >= 1 & msbb.raw$bbscore >= 4 & msbb.raw$NP.1 >= 2] = 'AD'

msbb <- msbb.raw %>%
  filter(RACE == 'W') %>%
  mutate(id = paste0('MSBB', WGS),
        study = 'MSBB',
        dx = recode(dx.raw, 'AD' = 'AD', 'CONTROL' = 'CTRL', 'OTHER' = NA_character_),
        dx = fct_relevel(dx, 'CTRL', 'AD'),
        apoe = recode(apoe, 'NANA' = NA_character_),
        apoe4 = recode(apoe, '22' = 'e4-', '23' = 'e4-', '33' = 'e4-', '24' = 'e4+', '34' = 'e4+', '44' = 'e4+'),
        aod = str_replace(AOD, '\\+', ''),
        aod = as.numeric(aod),
        aod_cat = cut(aod, c(50, 60, 70, 80, 90, Inf), c('50-69', '60-69', '70-79', '80-89', '90+'), right = FALSE),
        SEX = as.factor(SEX), dx = as.factor(dx), mtcn_avg = as.numeric(mtcn_avg), dx.raw = as.factor(dx.raw), apoe = as.factor(apoe), apoe4 = as.factor(apoe4)) %>%
  select(id, study, sex = SEX, apoe, apoe4, aod, aod_cat, mtcn_avg, dx.raw, dx)

## Mayo Clinic
mayo.raw <- read_tsv('data/AMPAD/mayo/WGS_Metadata.txt') %>%
  mutate(WGS_Participant_ID = as.character(WGS_Participant_ID)) %>%
  left_join(mtdna, by = c('WGS_Participant_ID' = 'sample'))

mayo <- mayo.raw %>%
  mutate(id = paste0('MAYO', WGS_Participant_ID),
        study = 'MAYO',
        apoe4 = recode(ApoE, '22' = 'e4-', '23' = 'e4-', '33' = 'e4-', '24' = 'e4+', '34' = 'e4+', '44' = 'e4+'),
        dx = recode(Diagnosis, 'AD' = 'AD', 'Control' = 'CTRL', 'Pathologic Aging' = NA_character_, 'PSP' = NA_character_),
        dx = fct_relevel(dx, 'CTRL', 'AD'),
        aod = str_replace(AgeAtDeath, '90_or_above', '90'),
        aod = as.numeric(aod),
        aod_cat = cut(aod, c(50, 60, 70, 80, 90, Inf), c('50-69', '60-69', '70-79', '80-89', '90+'), right = FALSE),
        Sex = as.factor(Sex), dx = as.factor(dx), mtcn_avg = as.numeric(mtcn_avg), Diagnosis = as.factor(Diagnosis), ApoE = as.factor(ApoE), apoe4 = as.factor(apoe4)) %>%
  select(id, study, sex = Sex, apoe = ApoE, apoe4, aod, aod_cat, mtcn_avg, dx.raw = Diagnosis, dx)


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

  rosmap %>%
    mutate(dx = recode(dx, '1' = 'CTRL', '4' = 'AD', '5' = 'AD'))

  recode(rosmap$dx, '1' = 'CTRL', '2' = NA)

  glm(dx ~ mtcn_avg + aod + sex + apoe4, data = rosmap, family = 'binomial') %>%
    tidy()

  glm(dx ~ mtcn_avg + aod + sex + apoe4, data = mayo, family = 'binomial') %>%
    summary()

  glm(dx ~ mtcn_avg + aod + sex + apoe4, data = msbb, family = 'binomial') %>%
    tidy()

  ampad <- rosmap %>%
    bind_rows(mayo) %>%
    bind_rows(msbb) %>%
    mutate(study = as_factor(study))

  glm(dx ~ mtcn_avg + aod + sex  + study, data = ampad, family = 'binomial') %>%
    tidy()
    glm(dx ~ mtcn_avg + aod + study, data = ampad, family = 'binomial') %>%
      tidy()
  t.test(mtcn_avg ~ dx, data = ampad)

  lm(mtcn_avg ~  dx + aod + sex + study, data = ampad)
  lm(aod ~  dx + mtcn_avg + sex + study, data = ampad) %>% tidy()
  m <- polr(aod_cat ~  dx + mtcn_avg + sex + study, data = ampad)
  (ctable <- coef(summary(m)))
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
  (ctable <- cbind(ctable, "p value" = p))
