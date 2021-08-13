brary(dplyr)

rosmap <- read_tsv('/sc/orga/projects/LOAD/Data/AMPAD/rosmap/WGS_Metadata.txt') %>% 
  select(bam_file) 
msbb <- read_tsv('/sc/orga/projects/LOAD/Data/AMPAD/msbb/WGS_Metadata.txt')  %>% 
  select(files.x) %>% rename(bam_file = files.x)
mayo <- read_tsv('/sc/orga/projects/LOAD/Data/AMPAD/mayo/WGS_Metadata.txt') %>% 
  select(files.x)  %>% rename(bam_file = files.x)

files <- rosmap %>% 
  bind_rows(msbb) %>% 
  bind_rows(mayo) %>% 
  mutate(bai = bam_file) %>% 
  mutate(bai = str_replace(bai, 'bam', 'bai')) %>% 
  gather('type', 'path')

sample_file <- str_split(files$path, pattern = '/', simplify = TRUE)[,10]
files$sample <- str_split(sample_file, pattern = '\\.', simplify = TRUE)[,1]

files %>% select(path) %>% write_tsv('/sc/orga/projects/LOAD/shea/Projects/mDNACN/data/bam_files.txt', col_names = F)
files %>% 
  filter(type == 'bam_file') %>%
  select(sample) %>% 
  write_tsv('/sc/orga/projects/LOAD/shea/Projects/mDNACN/data/sampleIDs.txt', col_names = F)




