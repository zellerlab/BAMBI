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
fn.meta <- here('data', 'meta_Zeevi_2014.tsv')
meta.zeevi <- read_tsv(fn.meta) %>%
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
meta.xie <- read_tsv(fn.meta) %>%
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
meta.goodrich <- read_tsv(fn.meta) %>%
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
meta.schirmer <- read_tsv(fn.meta) %>%
  filter(!is.na(Sample_ID)) %>%
  as.data.frame()
rownames(meta.schirmer) <- meta.schirmer$Sample_ID
meta.schirmer$Individual_ID <- meta.schirmer$Sample_ID
meta.schirmer$Timepoint <- 0

# ##############################################################################
# Milieu interieur
fn.feat <- here('data', 'motus_MilieuInterieur_WGS.tsv')
feat.mi.wgs <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                       check.names = FALSE, quote = '', comment.char = '',
                       row.names = 1, header = TRUE)
fn.meta <- here('data', 'meta_MilieuInterieur_WGS.tsv')
meta.mi.wgs <- read_tsv(fn.meta) %>%
  filter(!is.na(Sample_ID)) %>%
  as.data.frame()
rownames(meta.mi.wgs) <- meta.mi.wgs$Sample_ID
# 16S
fn.feat <- here('data', 'otus_MilieuInterieur_16S.tsv')
feat.mi.16S <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                          check.names = FALSE, quote = '', comment.char = '',
                          row.names = 1, header = TRUE)
fn.meta <- here('data', 'meta_MilieuInterieur_16S.tsv')
meta.mi.16s <- read_tsv(fn.meta) %>%
  filter(!is.na(Sample_ID)) %>%
  as.data.frame()
rownames(meta.mi.16s) <- meta.mi.16s$Sample_ID


dataset.list <- list(
  'Zeevi_WGS'=list('feat'=feat.zeevi, 'meta'=meta.zeevi),
  'Schirmer_WGS'=list('feat'=feat.schirmer, 'meta'=meta.schirmer),
  'TwinsUK_WGS'=list('feat'=feat.xie, 'meta'=meta.xie),
  'TwinsUK_16S'=list('feat'=feat.goodrich, 'meta'=meta.goodrich),
  'MI_WGS'=list('feat'=feat.mi.wgs, 'meta'=meta.mi.wgs),
  'MI_16S'=list('feat'=feat.mi.16S, 'meta'=meta.mi.16s),
  'Zeevi_KEGG'=list('feat'=feat.zeevi.kegg, 'meta'=meta.zeevi))
