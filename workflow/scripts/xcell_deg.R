## ========================================================================== ##
## Evaluate if gene signatures used by xCell for Neurons and Astrocytes are 
## Differentialy expressed in AMP-AD AD cases and controls
## DEG may bias xCell scores for Neurons and Astrocytes
## ========================================================================== ##

library(synapser)
library(tidyverse)

synLogin()

# xCell Gene signatures
## https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-017-1349-1/MediaObjects/13059_2017_1349_MOESM3_ESM.xlsx
xcell <- readxl::read_xlsx("raw/xCell/GeneSignatures.xlsx") %>% 
  select(-`# of genes`) %>% 
  pivot_longer(., -Celltype_Source_ID, names_to = "col", values_to = "gene") %>% 
  filter(!is.na(gene)) %>%
  filter(str_detect(Celltype_Source_ID, "^Neuro|^Astro")) %>%
  mutate(col = case_when(str_detect(Celltype_Source_ID, "^Astro") ~ "Astrocytes", 
                         str_detect(Celltype_Source_ID, "^Neuro") ~ "Neurons")) %>% 
  group_by(col) %>%
  distinct(gene)

## Number of Gene Signatures for Neurons and Astrocytes
xcell %>% count(col)  

# AMP-AD DEG 
## https://www.cell.com/cell-reports/fulltext/S2211-1247(20)30889-5
deg_syn <- synGet("syn14237651") 
de <- read_tsv(deg_syn$path) 

xcell_de <- left_join(xcell, de, by = c("gene" = "hgnc_symbol"))

## Restrict to AD vs control; all sexs; harmonized diagnosis; ROSMAP
xcell_de %>% 
  ungroup() %>%
  filter(Comparison == "AD-CONTROL") %>% 
  filter(Sex == "ALL") %>%
  filter(Model == "Diagnosis") %>%
  filter(Study == "ROSMAP") %>%
  relocate(Study, before = Model) %>%
  group_by(col) %>%
  filter(adj.P.Val < 0.05) %>% 
  select(Study, Model, col, gene, Tissue, logFC, CI.L, CI.R, AveExpr, Direction, P.Value, adj.P.Val) %>%
#  distinct(gene) %>%
  count(col, Direction) %>% 
  print(n = Inf)
