# ##############################################################################
#
## Load the data
#
# ##############################################################################

library("here")
library("tidyverse")

# ##############################################################################
# Zeevi
fn.feat <- here('data', 'motus_Zeevi_2014.tsv')
feat.zeevi <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                         check.names = FALSE, quote = '', comment.char = '',
                         row.names = 1, header = TRUE)
fn.feat <- here('data', 'KEGG_Zeevi_2014.tsv')
feat.zeevi.kegg <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                              check.names = FALSE, quote = '', 
                              comment.char = '',
                              row.names = 1, header = TRUE)
# transform KEGG to integer values
feat.zeevi.kegg <- round(feat.zeevi.kegg)
feat.zeevi.kegg[] <- lapply(feat.zeevi.kegg, as.integer)
feat.zeevi.kegg <- feat.zeevi.kegg[-which(rownames(feat.zeevi.kegg) == '-1'),]
feat.zeevi.kegg <- feat.zeevi.kegg[,-which(colSums(feat.zeevi.kegg)< 250000)]
fn.meta <- here('data', 'meta_Zeevi_2014.tsv')
meta.zeevi <- read_tsv(fn.meta, col_types = cols()) %>%
  filter(!is.na(Sample_ID)) %>%
  as.data.frame()
rownames(meta.zeevi) <- meta.zeevi$Sample_ID

# ##############################################################################
# TwinsUK
fn.feat <- here('data', 'motus_Xie_2016.tsv')
feat.xie <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                       check.names = FALSE, quote = '', comment.char = '',
                       row.names = 1, header = TRUE)
fn.meta <- here('data', 'meta_Xie_2016.tsv')
meta.xie <- read_tsv(fn.meta, col_types = cols()) %>%
  filter(!is.na(Sample_ID)) %>%
  as.data.frame()
rownames(meta.xie) <- meta.xie$Sample_ID
meta.xie$Individual_ID <- meta.xie$Sample_ID
meta.xie$Timepoint <- 0
# 16S
fn.feat <- here('data', 'otus_Goodrich_2014.tsv')
feat.goodrich <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                            check.names = FALSE, quote = '', comment.char = '',
                            row.names = 1, header = TRUE)
fn.meta <- here('data', 'meta_Goodrich_2014.tsv')
meta.goodrich <- read_tsv(fn.meta, col_types = cols()) %>%
  filter(!is.na(Sample_ID)) %>%
  filter(DiseaseState != 'OB') %>% 
  select(Sample_ID, age, body_mass_index) %>% 
  as.data.frame()
rownames(meta.goodrich) <- meta.goodrich$Sample_ID
meta.goodrich$Individual_ID <- meta.goodrich$Sample_ID
meta.goodrich$Timepoint <- 0

# ##############################################################################
# Schirmer
fn.feat <- here('data', 'motus_Schirmer_2016.tsv')
feat.schirmer <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                            check.names = FALSE, quote = '', comment.char = '',
                            row.names = 1, header = TRUE)
fn.meta <- here('data', 'meta_Schirmer_2016.tsv')
meta.schirmer <- read_tsv(fn.meta, col_types = cols()) %>%
  filter(!is.na(Sample_ID)) %>%
  as.data.frame()
rownames(meta.schirmer) <- meta.schirmer$Sample_ID
meta.schirmer$Individual_ID <- meta.schirmer$Sample_ID
meta.schirmer$Timepoint <- 0


# ##############################################################################
# HMP
fn.meta <- list.files(here('data'), pattern='meta_HMP_', full.names = TRUE)
meta.all <- map(fn.meta, read_tsv)
names(meta.all) <- str_remove(fn.meta, '.*HMP_') %>% str_remove('.tsv')
fn.feat <- list.files(here('data'), pattern='otus_HMP_', full.names = TRUE)
feat.all <- map(fn.meta, .f = function(x){
  read.table(x, sep='\t', stringsAsFactors = FALSE, quote = '', 
             comment.char = '', header = TRUE, row.names = 1)})
names(feat.all) <- str_remove(fn.feat, '.*HMP_') %>% str_remove('.tsv')


dataset.list <- list(
  'Zeevi_WGS'=list('feat'=feat.zeevi, 'meta'=meta.zeevi),
  'Schirmer_WGS'=list('feat'=feat.schirmer, 'meta'=meta.schirmer),
  'TwinsUK_WGS'=list('feat'=feat.xie, 'meta'=meta.xie),
  'TwinsUK_16S'=list('feat'=feat.goodrich, 'meta'=meta.goodrich),
  'Zeevi_KEGG'=list('feat'=feat.zeevi.kegg, 'meta'=meta.zeevi),
  'HMP_airways'=list('feat'=feat.all$airways, 'meta'=meta.all$airways),
  'HMP_saliva'=list('feat'=feat.all$saliva, 'meta'=meta.all$saliva),
  'HMP_skin'=list('feat'=feat.all$skin, 'meta'=meta.all$skin),
  'HMP_stool'=list('feat'=feat.all$stool, 'meta'=meta.all$stool),
  'HMP_vagina'=list('feat'=feat.all$vagina, 'meta'=meta.all$vagina))
