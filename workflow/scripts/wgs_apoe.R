library(tidyverse)
library(vcfR)

## ROSMAP
system("bcftools view -r 19:45162079-45662079 -Ov -o DEJ_11898_B01_GRM_WGS_2017-05-15_APOE.recalibrated_variants.vcf.gz /sc/arion/projects/AMPADWGS/RawData/DEJ11898_06_16_2017/Project_DEJ_11898_B01_GRM_WGS.JGvariants.2017-06-06/jgwd/joint_vcf/DEJ_11898_B01_GRM_WGS_2017-05-15_19.recalibrated_variants.vcf.gz")

rosmap.wgsqc <- read_csv("data/AMPAD_extra/rosmap/WGS_sample_QC_info.csv", guess_max = 10000)
vcf <- vcfR::read.vcfR("data/AMPAD_extra/rosmap/DEJ_11898_B01_GRM_WGS_2017-05-15_APOE.recalibrated_variants.vcf.gz")

Z <- vcfR2tidy(vcf, format_fields = c("GT"))
apoe <- Z$gt %>%
  filter(POS %in% c(45412079, 45411941)) %>%
  mutate(snp = as_factor(POS),
         snp = fct_recode(snp, "rs7412" = "45412079", "rs429358" = "45411941")) %>%
  select(Indiv, snp, gt_GT_alleles) %>%
  pivot_wider(values_from = gt_GT_alleles, names_from = snp) %>%
    mutate(
      apoe = case_when(
          rs429358 == "T/C" & rs7412 == "C/T" ~ 24,
          rs429358 == "T/T" & rs7412 == "T/T" ~ 22,
          rs429358 == "T/T" & rs7412 == "C/T" ~ 23,
          rs429358 == "T/T" & rs7412 == "C/C" ~ 33,
          rs429358 == "T/C" & rs7412 == "C/C" ~ 34,
          rs429358 == "C/C" & rs7412 == "C/C" ~ 44
        )
    ) %>%
  left_join(select(rosmap.wgsqc, projid, WGS_id), by = c('Indiv' = 'WGS_id')) %>%
  distinct(projid, .keep_all = TRUE) %>%
  select(-Indiv)

write_tsv(apoe, "data/AMPAD_extra/rosmap/wgs_apoe.tsv")

## mayo
system("bcftools view -r 19:45162079-45662079 -Ov -o NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_APOE.recalibrated_variants.Mayo.vcf.gz NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_19.recalibrated_variants.Mayo.vcf.gz")

vcf <- vcfR::read.vcfR(
  "data/AMPAD_extra/mayo/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_APOE.recalibrated_variants.Mayo.vcf.gz"
  )

Z <- vcfR::vcfR2tidy(vcf, format_fields = c("GT"))
apoe <- Z$gt %>%
  filter(POS %in% c(45412079, 45411941)) %>%
  mutate(snp = as_factor(POS),
         snp = fct_recode(snp, "rs7412" = "45412079", "rs429358" = "45411941")) %>%
  select(Indiv, snp, gt_GT_alleles) %>%
  pivot_wider(values_from = gt_GT_alleles, names_from = snp) %>%
    mutate(
      apoe = case_when(
          rs429358 == "T/C" & rs7412 == "C/T" ~ 24,
          rs429358 == "T/T" & rs7412 == "T/T" ~ 22,
          rs429358 == "T/T" & rs7412 == "C/T" ~ 23,
          rs429358 == "T/T" & rs7412 == "C/C" ~ 33,
          rs429358 == "T/C" & rs7412 == "C/C" ~ 34,
          rs429358 == "C/C" & rs7412 == "C/C" ~ 44
        )
    )

write_tsv(apoe, "data/AMPAD_extra/mayo/wgs_apoe.tsv")
