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
# CRC and IBD (take only the CD subtype) datasets

# CRC
crc.meta <- list()
crc.feat <- list()
for (i in c('Zeller_2014', 'Feng_2015', 'Yu_2017', 
            'Vogtmann_2016', 'Wirbel_2019')){
  meta <- read_tsv(paste0(data.location, 'metadata/meta_', i, '.tsv'),
                   col_types=cols()) %>% 
    select(Sample_ID, Group, Sampling_rel_to_colonoscopy) %>% 
    mutate(Study=i) %>% 
    filter(Group %in% c('CTR', 'CRC'))
  
  feat <- read.table(paste0(data.location, 
                            'tax_profiles/mOTUs_v2.5/', i, '.motus'),
                     sep='\t', stringsAsFactors = FALSE,
                     check.names = FALSE, quote = '', comment.char = '',
                     row.names = 1, header = TRUE, skip=61)
  stopifnot(all(meta$Sample_ID %in% colnames(feat)))
  feat <- feat[,colSums(feat) > 100]
  meta <- meta %>% filter(Sample_ID %in% colnames(feat))
  
  crc.feat[[length(crc.feat)+1]] <- feat[,meta$Sample_ID]
  crc.meta[[i]] <- meta
  message("Finished dataset ", i)
}
crc.meta <- bind_rows(crc.meta)
feat <- do.call(cbind, crc.feat)
feat.rel <- prop.table(as.matrix(feat), 2)

# prevalence filter
prev.mat <- vapply(unique(crc.meta$Study), FUN = function(s){
    rowMeans(feat[,crc.meta %>% filter(Study==s) %>% pull(Sample_ID)] != 0)}, 
    FUN.VALUE = double(nrow(feat.rel)))
f.idx <- names(which(rowSums(prev.mat > 0.05) > 2))
feat.rel <- feat.rel[f.idx,]

# save tables
write_tsv(crc.meta, here('data', 'meta_crc.tsv'))
write.table(feat, here('data', 'motus_crc_meta.tsv'),
            sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(feat.rel, here('data', 'motus_crc_rel_meta.tsv'),
            sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)


# CD
cd.meta <- list()
cd.feat <- list()
for (i in c('Lewis_2015', 'He_2017', 'metaHIT', 'HMP2', 'Franzosa_2019')){
  meta <- read_tsv(paste0(temp.loc, '/metadata/meta_', i, '.tsv'),
                   col_types=cols())
  if (i == 'Lewis_2015'){
    meta <- meta %>% 
      group_by(Subject) %>% 
      filter(Time==min(Time)) %>% 
      ungroup()
  } else if (i == 'metaHIT'){
    meta <- meta %>% 
      filter(Country=='spanish') %>% 
      group_by(Individual_ID) %>% 
      filter(Sampling_day==min(Sampling_day)) %>% 
      ungroup()
  } else if (i == 'HMP2'){
    meta <- meta %>% 
      mutate(Group=case_when(Group=='nonIBD'~'CTR', TRUE~Group)) %>% 
      group_by(Individual_ID) %>% 
      filter(Timepoint==min(Timepoint)) %>% 
      ungroup()
  }
  
  meta <- meta %>% select(Sample_ID, Group) %>% 
    mutate(Study=i) %>% 
    mutate(Sample_ID=make.names(Sample_ID)) %>% 
    filter(Group %in% c('CTR', 'CD'))
  
  feat <- read.table(paste0(temp.loc, 
                            '/motus/', i, '_motus.tsv'),
                     sep='\t', stringsAsFactors = FALSE,
                     check.names = FALSE, quote = '', comment.char = '',
                     row.names = 1, header = TRUE)
  colnames(feat) <- make.names(colnames(feat))
  feat <- feat[,colSums(feat) > 100]
  meta <- meta %>% filter(Sample_ID %in% colnames(feat))

  cd.feat[[length(cd.feat)+1]] <- feat[,meta$Sample_ID]
  cd.meta[[i]] <- meta
  message("Finished dataset ", i)
}
cd.meta <- bind_rows(cd.meta)
feat <- do.call(cbind, cd.feat)
feat.rel <- prop.table(as.matrix(feat), 2)

# also genus level
tax.info <- read_tsv(here('data', 'motus_taxonomy.tsv'))
tax.info <- tax.info %>% 
  filter(!str_detect(family, '^NA'))
bins.unique <- unique(tax.info$genus)
feat.genus <- matrix(NA, nrow=length(bins.unique), ncol=ncol(feat),
                     dimnames=list(bins.unique, colnames(feat)))
pb <- progress_bar$new(total=length(bins.unique))
feat.temp <- feat
rownames(feat.temp) <- str_extract(rownames(feat), '((ref|meta)_mOTU_v25_[0-9]{5}|^-1$)')
for (x in bins.unique){
  feat.genus[x,] <- colSums(feat.temp[tax.info %>% 
                                      filter(genus==x) %>% 
                                      pull(mOTU_ID),,drop=FALSE])
  pb$tick()
}
feat.genus <- rbind(feat.genus, feat.temp['-1',])
# prevalence filter
prev.mat <- vapply(unique(cd.meta$Study), FUN = function(s){
  rowMeans(feat[,cd.meta %>% filter(Study==s) %>% pull(Sample_ID)] != 0)}, 
  FUN.VALUE = double(nrow(feat.rel)))
f.idx <- names(which(rowSums(prev.mat > 0.05) > 2))
feat.rel <- feat.rel[f.idx,]

prev.mat <- vapply(unique(cd.meta$Study), FUN = function(s){
  rowMeans(feat.genus[,cd.meta %>% filter(Study==s) %>% pull(Sample_ID)] != 0)}, 
  FUN.VALUE = double(nrow(feat.genus)))
f.idx <- names(which(rowSums(prev.mat > 0.05) > 2))
feat.genus <- feat.genus[f.idx,]

# save tables
write_tsv(cd.meta, here('data', 'meta_cd.tsv'))
write.table(feat, here('data', 'motus_cd_meta.tsv'),
            sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(feat.rel, here('data', 'motus_cd_rel_meta.tsv'),
            sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(feat.genus, here('data', 'motus_cd_genus.tsv'),
            sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
