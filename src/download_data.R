# ##############################################################################
#
## Download and prep data for the project
#
# ##############################################################################

renv::load(here::here())

library("tidyverse")
library("here")

data.location <- 'https://intranet.embl.de/download/zeller/'

# Zeevi data
fn.feat <- paste0(data.location, 'tax_profiles/mOTUs_v2.5/Zeevi.motus')
feat <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                   check.names = FALSE, quote = '', comment.char = '',
                   row.names = 1, header = TRUE, skip=61)
write.table(feat, here('data', 'motus_Zeevi_2014.tsv'), sep='\t', 
            quote = FALSE, col.names = TRUE, row.names = TRUE)
fn.meta <- paste0(data.location, 'metadata/meta_Zeevi.tsv')
meta <- read_tsv(fn.meta) %>% 
  mutate(Gender=case_when(Gender==0~'F', Gender==1~'M'))
write_tsv(meta, here('data', 'meta_Zeevi_2014.tsv'))

# TwinsUK data
fn.feat <- paste0(data.location, 'tax_profiles/mOTUs_v2.5/TwinsUK.motus')
feat <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                   check.names = FALSE, quote = '', comment.char = '',
                   row.names = 1, header = TRUE, skip=61)
write.table(feat, here('data', 'motus_Xie_2016.tsv'), sep='\t', 
            quote = FALSE, col.names = TRUE, row.names = TRUE)
fn.meta <- paste0(data.location, 'metadata/meta_TwinsUK.tsv')
meta <- read_tsv(fn.meta)
write_tsv(meta, here('data', 'meta_Xie_2016.tsv'))

# metaCARDIS 
# TODO

# AGP
fn.feat <- paste0(data.location, 'tax_profiles/otus_AGP.tsv')
feat <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                   check.names = FALSE, quote = '', comment.char = '',
                   row.names = 1, header = TRUE)
write.table(feat, here('data', 'otus_AGP_2020.tsv'), sep='\t', 
            quote = FALSE, col.names = TRUE, row.names = TRUE)
fn.meta <- paste0(data.location, 'metadata/meta_AGP.tsv')
meta <- read_tsv(fn.meta)
write_tsv(meta, here('data', 'meta_AGP_2020.tsv'))
