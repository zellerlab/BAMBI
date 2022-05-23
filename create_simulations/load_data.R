# ##############################################################################
#
## Load the data
#
# ##############################################################################

library("here")
library("tidyverse")

fn.feat <- here('data', 'motus_Zeevi_2014.tsv')
feat.zeevi <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                         check.names = FALSE, quote = '', comment.char = '',
                         row.names = 1, header = TRUE)
fn.meta <- here('data', 'meta_Zeevi_2014.tsv')
meta.zeevi <- read_tsv(fn.meta) %>%
  filter(!is.na(Sample_ID)) %>%
  as.data.frame()
rownames(meta.zeevi) <- meta.zeevi$Sample_ID
# TwinsUK
fn.feat <- here('data', 'motus_Xie_2016.tsv')
feat.twinsuk <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                           check.names = FALSE, quote = '', comment.char = '',
                           row.names = 1, header = TRUE)
fn.meta <- here('data', 'meta_Xie_2016.tsv')
meta.twinsuk <- read_tsv(fn.meta) %>%
  filter(!is.na(Sample_ID)) %>%
  as.data.frame()
rownames(meta.twinsuk) <- meta.twinsuk$Sample_ID
meta.twinsuk$Individual_ID <- meta.twinsuk$Sample_ID
meta.twinsuk$Timepoint <- 0
