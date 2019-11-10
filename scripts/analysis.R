library(readxl)
library(forcats)
library(tidyverse)
library(broom)

#setwd('/Users/sheaandrews/LOAD_minerva/dummy/shea/Projects/mDNACN')

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
rosmpa_rna <- read_tsv('data/AMPAD/rosmap/RNAseq_WGS_Metadata.txt')
rosmap.raw <- read_tsv('data/AMPAD/rosmap/WGS_Metadata.txt') %>%
  left_join(mtdna, by = c('wgs_id' = 'sample')) %>%
  left_join(select(rosmpa_rna, wgs_id, rnaseq_id), by = 'wgs_id')

rosmap <- rosmap.raw %>%
  filter(race == 1) %>%
  mutate(aod = str_replace(age_death, '\\+', ''),
        aod = round(as.numeric(aod), 0),
        dx = recode(cogdx, '1' = 'CTRL', '4' = 'AD', '5' = 'AD'),
        dx = fct_relevel(dx, 'CTRL', 'AD'),
        apoe4 = recode(apoe_genotype, '22' = 'e4-', '23' = 'e4-', '33' = 'e4-', '24' = 'e4+', '34' = 'e4+', '44' = 'e4+'),
        aod_cat = cut(aod, c(50, 60, 70, 80, 90, Inf), c('50-69', '60-69', '70-79', '80-89', '90+'), right = FALSE),
        msex = as.factor(msex), dx = as.factor(dx), mtcn_avg = as.numeric(mtcn_avg),
        cogdx = as.factor(cogdx), apoe_genotype = as.factor(apoe_genotype), apoe4 = as.factor(apoe4),
        z_mtdnacn = scale(mtcn_avg, center = TRUE, scale = TRUE)[,1]) %>%
  select(id = wgs_id, rnaseq_id, study, sex = msex, apoe = apoe_genotype, apoe4, aod, aod_cat, mtcn_avg, z_mtdnacn, dx.raw = cogdx, dx)

## Mount Sinai Brain Bank
msbb.raw <- read_tsv('data/AMPAD/msbb/WGS_Metadata.txt') %>%
#  unite(apoe, c('Apo1', 'Apo2'), sep = "") %>%
  mutate(WGS = as.character(WGS)) %>%
  left_join(mtdna, by = c('WGS' = 'sample'))

msbb.raw$dx.raw = 'OTHER'
msbb.raw$dx.raw[msbb.raw$CDR <= 0.5 & msbb.raw$bbscore <= 3 & msbb.raw$NP.1 <= 1] = 'CONTROL'
msbb.raw$dx.raw[msbb.raw$CDR >= 1 & msbb.raw$bbscore >= 4 & msbb.raw$NP.1 >= 2] = 'AD'

apoe <- read_table2('apoe_only.recode.vcf', comment = '##') %>%
  rename(CHROM = `#CHROM`) %>%
  filter(POS %in% c(45412079, 45411941)) %>%
  mutate(ID = paste(CHROM, POS, REF, ALT, sep = ":")) %>%
  select(-CHROM, -POS, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>%
  gather('sample', 'genotype', 2:ncol(.)) %>%
  separate(genotype, c('A1', 'A2')) %>%
  mutate(A1 = as.numeric(A1),
         A2 = as.numeric(A2),
         rsid = ifelse(ID == '19:45411941:T:C', 'rs429358', 'rs7412')) %>%
  rowwise() %>%
  mutate(alleles = sum(A1, A2),
         sample = paste0('MSBB', sample)) %>%
  select(-A1, -A2, -ID) %>%
  ## 19:45411941:T:C, 19:45412079:C:T
  spread(rsid, alleles) %>%
  mutate(A1 = recode(rs429358, '0' = '3', '1' = '4', '2' = '4'),
         A2 = recode(rs7412, '0' = '3', '1' = '2', '2' = '2')) %>%
  unite(apoe, c('A2', 'A1'), sep = "")

msbb <- msbb.raw %>%
  filter(RACE == 'W') %>%
  mutate(id = paste0('MSBB', WGS)) %>%
  left_join(apoe, by = c('id' = 'sample')) %>%
  mutate(study = 'MSBB',
        dx = recode(dx.raw, 'AD' = 'AD', 'CONTROL' = 'CTRL', 'OTHER' = NA_character_),
        dx = fct_relevel(dx, 'CTRL', 'AD'),
        cdr.dx = recode(CDR, '0' = '0', '0.5' = '0', '1' = '1', '2' = '1', '3' = '1', '4' = '1', '5' = '1'),
        apoe4 = recode(apoe, '22' = 'e4-', '23' = 'e4-', '33' = 'e4-', '24' = 'e4+', '34' = 'e4+', '44' = 'e4+'),
        aod = str_replace(AOD, '\\+', ''),
        aod = as.numeric(aod),
        aod_cat = cut(aod, c(50, 60, 70, 80, 90, Inf), c('50-69', '60-69', '70-79', '80-89', '90+'), right = FALSE),
        SEX = as.factor(SEX), dx = as.factor(dx), mtcn_avg = as.numeric(mtcn_avg),
        dx.raw = as.factor(dx.raw), cdr.dx = as.factor(cdr.dx),
        apoe = as.factor(apoe), apoe4 = as.factor(apoe4),
        z_mtdnacn = scale(mtcn_avg, center = TRUE, scale = TRUE)[,1]) %>%
  select(id, study, sex = SEX, apoe, apoe4, aod, aod_cat, mtcn_avg, z_mtdnacn, dx.raw, cdr.dx, dx)

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
        Sex = as.factor(Sex), dx = as.factor(dx), mtcn_avg = as.numeric(mtcn_avg),
        Diagnosis = as.factor(Diagnosis), ApoE = as.factor(ApoE), apoe4 = as.factor(apoe4),
        z_mtdnacn = scale(mtcn_avg, center = TRUE, scale = TRUE)) %>%
  select(id, study, sex = Sex, apoe = ApoE, apoe4, aod, aod_cat, mtcn_avg, z_mtdnacn, dx.raw = Diagnosis, dx)


glm(dx ~ mtcn_avg + aod + sex + apoe4, data = rosmap, family = 'binomial') %>%
 tidy()

glm(dx ~ mtcn_avg + aod + sex + apoe4, data = mayo, family = 'binomial') %>%
  summary()

glm(dx ~ mtcn_avg + aod + sex + apoe4, data = msbb, family = 'binomial') %>%
  tidy()
res <- glm(cdr.dx ~ z_mtdnacn + aod + sex + apoe4, data = msbb, family = 'binomial')
res %>% tidy()
res %>% tidy()


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

write_tsv(output/)
## gplot
ggplot(data = msbb, aes(x = cdr.dx, y = mtcn_avg, colour = cdr.dx, fill = cdr.dx)) +
  geom_quasirandom(width = 0.3, alpha = 0.5, size = 0.5) + #facet_wrap(. ~ is.na(apoe4)) +
  geom_boxplot(width = 0.25, colour = 'black', alpha = 0.5, size = 0.5) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=8),
        axis.title.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  labs(y = 'mtDNA-CN') +
  scale_x_discrete(labels = c('Control', 'Case')) +
  scale_colour_manual(values = c('#0571b0', '#ca0020')) +
  scale_fill_manual(values = c('#0571b0', '#ca0020'))
ggsave('~/Dropbox/Research/Grants/AARF - 2019/Q4/Drafts/Figure1.tiff', units = 'cm', width = 5, height = 5)


ggplot(data = msbb, aes(x = as.factor(CDR), y = mtcn_avg, colour = cdr.dx)) +
  geom_quasirandom(width = 0.3, alpha = 0.5, size = 0.5) + #facet_wrap(. ~ is.na(apoe4)) +
  geom_boxplot(width = 0.25, colour = 'black', alpha = 0.5, size = 0.25, outlier.shape = NA) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        text = element_text(size=8),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(y = 'mtDNA-CN', x = 'CDR Score') +
  scale_colour_manual(values = c('#0571b0', '#ca0020'),
                      labels = c('Control', 'Case')) +
  scale_fill_manual(values = c('#0571b0', '#ca0020'))
ggsave('~/Dropbox/Research/Grants/AARF - 2019/Q4/Drafts/Figure1_v1.jpg', units = 'cm', width = 7.5, height = 5)
