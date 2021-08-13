library("reader")
setwd('~/LOAD_minerva/dummy/shea/Projects/mDNACN')

rosmap.raw <- read_tsv('data/AMPAD/rosmap/WGS_Metadata.txt') %>% 
  select(ID = wgs_id, study) %>% 
  mutate(ID = as.character(ID))
msbb.raw <- read_tsv('data/AMPAD/msbb/WGS_Metadata.txt') %>% 
  mutate(study = 'MSBB') %>% 
  select(ID = WGS, study) %>% 
  mutate(ID = as.character(ID))
mayo.raw <- read_tsv('data/AMPAD/mayo/WGS_Metadata.txt')  %>% 
  mutate(study = 'mayo') %>% 
  select(ID = WGS_Participant_ID, study) %>% 
  mutate(ID = as.character(ID))

ids <- bind_rows(rosmap.raw, msbb.raw, mayo.raw)

## ===================================================== ## 
##                      Haplogrep2
## ===================================================== ## 

haplogrep.files <- list.files(pattern = '.haplogrep.txt', recursive = T)
haplogrep <- haplogrep.files %>% 
  map(., function(x){
    message('\n', x, '\n')
    read_tsv(x, col_types = list(SampleID = col_character())) 
  }) %>% 
  bind_rows()
haplogrep$macro <- with(haplogrep,
                        ifelse(str_detect(Haplogroup, "^L"),
                               substr(Haplogroup, start = 1, stop = 2),
                               substr(Haplogroup, start = 1, stop = 1)
                        ))
haplogrep <- haplogrep %>% mutate(model = 'individual')

ggplot(haplogrep, aes(x = Quality, fill = macro)) + geom_histogram() + theme_bw() + xlim(0.4, 1.05) + ylim(0, 1000)


hj <- read_tsv('data/freebayes/joint/AMPAD.chrM.haplogrep.txt') %>% mutate(model = 'merge')
hj$macro <- with(haplogrep,
                        ifelse(str_detect(Haplogroup, "^L"),
                               substr(Haplogroup, start = 1, stop = 2),
                               substr(Haplogroup, start = 1, stop = 1)
                        ))

bind_rows(haplogrep, hj) %>% 
  ggplot(., aes(x = Quality, fill = macro)) + geom_histogram() + theme_bw() + facet_grid(rows = 'model')


test <- left_join(haplogrep, hj, by = 'SampleID')

count(test, Haplogroup.x == Haplogroup.y)  
filter(test, Haplogroup.x != Haplogroup.y)



## ===================================================== ## 
##                      PhyMer
## ===================================================== ## 

phymer.files <- list.files(pattern = '.phymer.txt', recursive = T)
hg <- phymer.files %>% 
  map(., function(x){
    message('\n', x, '\n')
    out <- read_csv(x, skip = 12, col_names = F) %>% 
      mutate(ID = str_extract(x, "(?<=phymer/).*(?=\\.phymer)"), 
             X1 = str_remove_all(X1, '[:punct:]'), 
             X4 = as.numeric(str_remove_all(X4, ']')))
    select(out, ID, hg = X1, X2, X3, X4)
  }) %>% 
  bind_rows() %>% 
  left_join(ids)
hg$macro <- with(hg,
                 ifelse(str_detect(hg, "^L"),
                        substr(hg, start = 1, stop = 2),
                        substr(hg, start = 1, stop = 1)
                 ))

test <- group_by(hg, ID) %>% arrange(ID, X4) %>% slice(1) %>% ungroup()

ggplot(test, aes(x = X4, fill = macro)) + geom_histogram() + theme_bw() + xlim(0.5, 1)


hg %>% group_by(ID) %>% distinct(macro) %>% ungroup() %>% count(ID) %>% arrange(n)


test %>% select(X4, macro, study) %>% mutate(macro = as.factor(macro)) %>% group_by(study) %>% skim()
test %>% select(X4, macro, study) %>% mutate(macro = as.factor(macro)) %>% group_by(study) %>% count(macro) %>% pivot_wider(names_from = study, values_from = n)

test %>% su

count(test, hg)
count(test, macro)
ggplot(test, aes(x = X4, fill = macro)) + geom_histogram() + xlim(0,1)

n.readLines('149.phymer.txt', skip = 12, n = 5) %>% 
  tibble::enframe() %>% 
  seperate(value, )
test <- read_tsv('MAP15387421.phymer.txt', skip = 12, col_names = F)
test <- read_csv('data/phymer/1000.phymer.txt', skip = 12, col_names = F) %>% 
  mutate(#ID = str_extract(x, "(?<=phymer/).*(?=\\.phymer)"), 
         X1 = str_remove_all(X1, '[:punct:]'), 
         X4 = as.numeric(str_remove_all(X4, ']')))

ad.genes1mb %>% filter(loci == '4p16.1') %>% distinct(gene_name)



countLines('data/phymer/1000.phymer.txt')

