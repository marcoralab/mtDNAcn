library(synapser)
library(xCell)
library(tidyverse)

synLogin()

## Synapse IDs for metadata 
msbb_meta_syn <- synGet("syn6100548") 
rosmap_meta_syn <- synGet("syn4300313")
rosmap_ids_syn <- synGet("syn3382527")
mayocbe_meta_syn <- synGet("syn5223705")  
mayotcx_meta_syn <- synGet("syn3817650")

## Wrangle ID information 
msbb_meta <- read_csv(msbb_meta_syn$path) %>% 
  filter(fileType == "bam") %>% 
  mutate(study = "MSBB") %>%
  select(SampleID = sampleIdentifier, study, Tissue = BrodmannArea, ID = individualIdentifier, batch, RIN)
rosmap_ids <- read_csv(rosmap_ids_syn$path) %>% 
  filter(!is.na(rnaseq_id)) %>% 
  select(ID = projid, SampleID = rnaseq_id) %>% 
  distinct(SampleID, .keep_all = TRUE)
rosmap_meta <- read_tsv(rosmap_meta_syn$path) %>% 
  select(SampleID = Sampleid, RIN = RINcontinuous, batch = Batch) %>% 
  left_join(rosmap_ids) %>%
  mutate(Tissue = "DLPFC", 
         batch = as.character(batch), 
         ID = as.character(ID), 
         study = "ROSMAP") %>% 
  group_by(SampleID) %>%
  top_n(1, batch) %>% 
  ungroup()
mayo_meta <- bind_rows(
    read_csv(mayocbe_meta_syn$path), 
    read_csv(mayotcx_meta_syn$path)) %>% 
  select(SampleID, RIN) %>% 
  mutate(study = "MAYO") %>%
  separate(SampleID, c("ID", "Tissue"), remove = FALSE) 

## =============================================================================
## Estimate TPM from raw count data 

## https://gist.github.com/slowkow/c6ab0348747f86e2748b
count2tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

## Synapse IDs for raw counts 
## https://www.synapse.org/#!Synapse:syn9702085
## https://www.synapse.org/#!Synapse:syn17010685 
msbb_counts_syn <- synGet("syn8691099") 
rosmap_counts_syn <- synGet("syn8691134")
mayocbe_counts_syn <- synGet("syn8690904")  
mayotcx_counts_syn <- synGet("syn8690799")


msbb_counts <- read_tsv(msbb_counts_syn$path)
rosmap_counts <- read_tsv(rosmap_counts_syn$path)
mayocbe_counts <- read_tsv(mayocbe_counts_syn$path)  
mayotcx_counts <- read_tsv(mayotcx_counts_syn$path)

## Join Counts data
counts <- msbb_counts %>% 
  full_join(rosmap_counts) %>% 
  full_join(mayocbe_counts) %>% 
  full_join(mayotcx_counts) %>%
  filter(str_detect(feature, "^ENSG")) %>% 
  separate(feature, c('ensembl_gene_id', 'position'), sep = '\\.') 

## Pull gene length and gene symbol from ensebble
ensembl = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
geneIDs <- biomaRt::getBM(attributes=c('ensembl_gene_id', 'chromosome_name', "start_position", "end_position",
                                       "hgnc_symbol", "external_gene_name", "external_gene_source"), 
                          filters = 'ensembl_gene_id', 
                          values = pull(counts, ensembl_gene_id), 
                          mart = ensembl) %>% 
  as_tibble() %>% 
  mutate(length = end_position - start_position) %>%
  arrange(ensembl_gene_id) 

counts_IDs <- geneIDs %>% 
  left_join(counts, by = "ensembl_gene_id") %>% 
  distinct(external_gene_name, .keep_all = TRUE) 

## Estimate TPM 
tpm <- counts_IDs %>% 
  select(-ensembl_gene_id, -position, -chromosome_name, -start_position, -end_position,
         -hgnc_symbol, -external_gene_name, -external_gene_source, -length) %>% 
  map_dfc(., ~count2tpm(.x, counts_IDs$length)) %>% 
  as.data.frame() %>%
  magrittr::set_rownames(counts_IDs$external_gene_name) 

## =============================================================================
## Estimate TPM by agregating transcript abundance (TPM)
##  https://www.biostars.org/p/329625/

## Synapse IDs for transcript abundance 
## https://www.synapse.org/#!Synapse:syn9702085
## https://www.synapse.org/#!Synapse:syn17010685 
msbb_txtpm_syn <- synGet("syn10507727") 
rosmap_txtpm_syn <- synGet("syn10507749")
mayocbe_txtpm_syn <- synGet("syn10499229")  
mayotcx_txtpm_syn <- synGet("syn10507725")


msbb_txtpm <- read_tsv(msbb_txtpm_syn$path)
rosmap_txtpm <- read_tsv(rosmap_txtpm_syn$path)
mayocbe_txtpm <- read_tsv(mayocbe_txtpm_syn$path)  
mayotcx_txtpm <- read_tsv(mayotcx_txtpm_syn$path)

txtpm <- msbb_txtpm %>% 
  full_join(rosmap_txtpm) %>% 
  full_join(mayocbe_txtpm) %>% 
  full_join(mayotcx_txtpm)  

tpm <- txtpm %>%
  mutate(gid = gsub("^[^\\|]+\\|([^\\|]+).+\\|([^\\]+)\\|[^\\]+\\|([^\\]+).+", "\\1|\\2", feature)) %>%
  group_by(gid) %>%
  summarise_at(vars(-feature), sum) %>%
  separate(gid, c('gid', 'gene'), sep = '\\|') %>% 
  separate(gid, c('gid', 'position'), sep = '\\.') %>% 
  select(-gid, -position) %>% 
  distinct(gene, .keep_all = TRUE) %>%
  as.data.frame() %>% 
  column_to_rownames("gene")

## =============================================================================
## Run xCell
xCell.res <- xCellAnalysis(tpm, rnaseq = TRUE, cell.types.use = c("Neurons", "Astrocytes"))
xCell.p <- xCellSignifcanceBetaDist(xCell.res) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "cells") %>% 
  as_tibble() %>% 
  magrittr::set_colnames(c("cells", colnames(xCell.res))) %>% 
  pivot_longer(., -cells, names_to = "SampleID", values_to = "xcell") %>% 
  pivot_wider(names_from = cells, values_from = xcell) %>% 
  rename(Neurons.pval = Neurons, Astrocytes.pval = Astrocytes)

## Wrgable xCell output and merge meta data 
xCell_tab <- xCell.res %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "cells") %>% 
  as_tibble() %>% 
  pivot_longer(., -cells, names_to = "SampleID", values_to = "xcell") %>% 
  pivot_wider(names_from = cells, values_from = xcell) %>% 
  left_join(xCell.p, by = "SampleID") %>%
  filter(SampleID != "feature") %>% 
  left_join(
    bind_rows(msbb_meta, rosmap_meta, mayo_meta), by = "SampleID")

write_csv(xCell_tab, "/sc/arion/projects/LOAD/shea/Projects/mtDNAcn/data/xcell/ampad_xCell.csv")






















