# ##############################################################################
#
## Download and prep data for the project
#
# ##############################################################################

library("tidyverse")
library("here")

data.location <- 'https://intranet.embl.de/download/zeller/'

# ##############################################################################
# Zeevi

# metadata
fn.meta <- paste0(data.location, 'metadata/meta_Zeevi.tsv')
meta <- read_tsv(fn.meta) %>% 
  mutate(Gender=case_when(Gender==0~'F', Gender==1~'M'))
write_tsv(meta, here('data', 'meta_Zeevi_2014.tsv'))

# tax profiles
fn.feat <- paste0(data.location, 'tax_profiles/mOTUs_v2.5/Zeevi.motus')
feat <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                   check.names = FALSE, quote = '', comment.char = '',
                   row.names = 1, header = TRUE, skip=61)
write.table(feat, here('data', 'motus_Zeevi_2014.tsv'), sep='\t', 
            quote = FALSE, col.names = TRUE, row.names = TRUE)

# functional profiles for Zeevi
fn.feat <- paste0(data.location, 'func_profiles/KEGG/KEGG_Zeevi.tsv')
feat <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                   check.names = FALSE, quote = '', comment.char = '',
                   row.names = 1, header = TRUE)
write.table(feat, here('data', 'KEGG_Zeevi_2014.tsv'), sep='\t', 
            quote = FALSE, col.names = TRUE, row.names = TRUE)

# ##############################################################################
# TwinsUK

# metadata
fn.meta <- paste0(data.location, 'metadata/meta_TwinsUK.tsv')
meta <- read_tsv(fn.meta)
write_tsv(meta, here('data', 'meta_Xie_2016.tsv'))

fn.meta <- paste0(data.location, 'metadata/meta_Goodrich_2014.tsv')
meta <- read_tsv(fn.meta)
write_tsv(meta, here('data', 'meta_Goodrich_2014.tsv'))

# tax profiles shotgun
fn.feat <- paste0(data.location, 'tax_profiles/mOTUs_v2.5/TwinsUK.motus')
feat <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                   check.names = FALSE, quote = '', comment.char = '',
                   row.names = 1, header = TRUE, skip=61)
write.table(feat, here('data', 'motus_Xie_2016.tsv'), sep='\t', 
            quote = FALSE, col.names = TRUE, row.names = TRUE)

# tax profiles 16S
fn.feat <- paste0(data.location, 'tax_profiles/16S/otus_Goodrich_2014.tsv')
feat <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                   check.names = FALSE, quote = '', comment.char = '',
                   row.names = 1, header = TRUE)
# combine at the genus level
taxa.sp.removed <- rownames(feat) %>% str_remove(';s__.*$')
all.genera <- taxa.sp.removed %>% unique()
new.feat <- matrix(0, nrow=length(all.genera), ncol=ncol(feat), 
                   dimnames = list(all.genera, colnames(feat)))
for (g in all.genera){
  new.feat[g,] <- colSums(feat[which(taxa.sp.removed == g),,drop=FALSE])
}
stopifnot(all(colSums(feat) == colSums(new.feat)))
write.table(new.feat, here('data', 'otus_Goodrich_2014.tsv'), sep='\t', 
            quote = FALSE, col.names = TRUE, row.names = TRUE)

# ##############################################################################
# Schirmer

# metadata
fn.meta <- paste0(data.location, 'metadata/meta_Schirmer.tsv')
meta <- read_tsv(fn.meta)
write_tsv(meta, here('data', 'meta_Schirmer_2016.tsv'))

# tax profiles
fn.feat <- paste0(data.location, 'tax_profiles/mOTUs_v2.5/Schirmer.motus')
feat <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                   check.names = FALSE, quote = '', comment.char = '',
                   row.names = 1, header = TRUE, skip=61)
write.table(feat, here('data', 'motus_Schirmer_2016.tsv'), sep='\t', 
            quote = FALSE, col.names = TRUE, row.names = TRUE)

# ##############################################################################
# Milieu Interieur

# metadata
fn.meta <- paste0(data.location, 'metadata/meta_MilieuInterieur_WGS.tsv')
meta <- read_tsv(fn.meta) %>% 
  mutate(Timepoint=str_remove(Timepoint, '^V')) %>% 
  mutate(Timepoint=as.double(Timepoint))
write_tsv(meta, here('data', 'meta_MilieuInterieur_WGS.tsv'))
fn.meta <- paste0(data.location, 'metadata/meta_MilieuInterieur_16S.tsv')
meta <- read_tsv(fn.meta) %>% 
  mutate(Timepoint=str_remove(Timepoint, '^V')) %>% 
  mutate(Timepoint=as.double(Timepoint))
write_tsv(meta, here('data', 'meta_MilieuInterieur_16S.tsv'))

# tax profiles shotgun
fn.feat <- paste0(data.location, 'tax_profiles/mOTUs_v2.5/Milieu_interieur.tsv')
feat <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                   check.names = FALSE, quote = '', comment.char = '',
                   row.names = 1, header = TRUE, skip=2)
write.table(feat, here('data', 'motus_MilieuInterieur_WGS.tsv'), sep='\t', 
            quote = FALSE, col.names = TRUE, row.names = TRUE)

# tax profiles 16S
fn.feat <- paste0(data.location, 'tax_profiles/16S/otus_Milieu_interieur.tsv')
feat <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                   check.names = FALSE, quote = '', comment.char = '',
                   row.names = 1, header = TRUE, skip=2)
colnames(feat) <- str_remove(colnames(feat), '.mpseq')
# add number of unresolved features at genus level
all.reads <- read.table(fn.feat, sep='\t', stringsAsFactors = FALSE,
                        check.names = FALSE, quote = '', comment.char = '',
                        row.names=1, header = FALSE, nrow=2, skip=1)
colnames(all.reads) <- all.reads[2,]
colnames(all.reads) <- str_remove(colnames(all.reads), '.mpseq')
stopifnot(all(colnames(feat) == colnames(all.reads)))
all.reads <- as.numeric(all.reads[1,])
names(all.reads) <- colnames(feat)
feat <- rbind(feat, all.reads - colSums(feat))
rownames(feat)[nrow(feat)] <- 'Unassigned'

write.table(feat, here('data', 'otus_MilieuInterieur_16S.tsv'), sep='\t', 
            quote = FALSE, col.names = TRUE, row.names = TRUE)

# ##############################################################################
# TODO
# any other 16S datasets?
