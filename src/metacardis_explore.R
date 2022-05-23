# ##############################################################################
#
## Exploratory data analysis on metacardis dataset(s)
#
# ##############################################################################
library(tidyverse)
library(here)

## METADATA ----
meta.vars <- here('data','metacardis-2021', 'hub.metadata.variables.v10.r') %>%
  read_tsv()
meta <- here('data','metacardis-2021', 'hub.metadata.samples.v10.r') %>%
  read_tsv() %>%
  mutate(LABEL = fct_collapse(PATGROUPFINAL_C,
                              Healthy = c('8','1 RH'),
                              MetSyn = '1',
                              Obese = '2a',
                              ObeseBariatric = '2b',
                              T2D = '3',
                              CAD = '4',
                              CCAD = '5',
                              CCAD_HF = '6',
                              HF = '7')) %>%
  select(SampleID, LABEL, ANTIBIOTIC, ANTIBIOTICS_TOTAL, BETA_BL_C, 
         METFORMIN_C, DIABSTATUS_C, STATINE_C, contains('PPI'),
         food_37_other_alcoholic_beverages) %>%
  rename(ALCOHOL_G = 'food_37_other_alcoholic_beverages')

## MULTIPLE input data types -- no raw counts with clean names exists ----
# marker gene length scaled normalized fraction; anonymized
var.names <-  here('data','metacardis-2021',
                   'hub.taxon.motus.ID.read.scaled.fraction.variables.v3.r') %>%
  read_tsv() %>% select(VariableID, DisplayName)
#2194 samples available, 2036 baseline
motus.feat <- here('data','metacardis-2021',
                   'hub.taxon.motus.ID.read.scaled.fraction.mgsamples.v3.r') %>%
  read_tsv() %>%
  filter(str_detect(MGSampleID, 'M0')) %>%
  gather('VariableID','val',-MGSampleID) %>%
  mutate(MGSampleID = str_remove_all(MGSampleID, 'M0_')) %>%
  filter(str_detect(MGSampleID, 
                    paste(meta$SampleID, collapse = '|'))) %>%
  left_join(var.names) %>%
  select(-VariableID) %>%
  arrange(MGSampleID) %>%
  spread('MGSampleID','val') %>%
  rename(feature = 'DisplayName')

# gene length normalized fraction; anonymized
# var.names <-  here('data','metacardis-2021',
#                    'hub.taxon.MGS.fpkm.variables.v5.r') %>%
#   read_tsv() %>% select(VariableID, DisplayName)
# mgs.feat <- here('data','metacardis-2021',
#                  'hub.taxon.MGS.fpkm.mgsamples.v5.r') %>% read_tsv() %>%
#   gather('VariableID','val',-MGSampleID) %>%
#   left_join(var.names) %>%
#   select(-VariableID) %>%
#   arrange(MGSampleID) %>%
#   spread('MGSampleID','val') %>%
#   rename(feature = 'DisplayName')

# downsampled counts
# var.names <-  here('data','metacardis-2021',
#                    'hub.taxon.motus.ID.down.15000000.unscaled.variables.v3.r') %>%
#   read_tsv() %>% select(VariableID, DisplayName)
# motus.ds.feat <- here('data','metacardis-2021',
#                       'hub.taxon.motus.ID.down.15000000.unscaled.mgsamples.v3.r') %>%
#   read_tsv() %>%
#   gather('VariableID','val',-MGSampleID) %>%
#   left_join(var.names) %>%
#   select(-VariableID) %>%
#   spread('MGSampleID','val') %>%
#   rename(feature = 'DisplayName')

# incorrectly named from hub -- these are TSS normalized !
# fine for our purposes -- we are not simulating with these
var.names <- here('data','metacardis-2021',
                  'hub.taxon.motus.Genus.raw.variables.v3.r') %>%
  read_tsv() %>% select(VariableID, DisplayName)
# 134 genera and 2028 individuals @baseline
motus.genus <- here('data','metacardis-2021',
                    'hub.taxon.motus.Genus.raw.mgsamples.v3.r') %>%
  read_tsv() %>%
  filter(str_detect(MGSampleID, 'M0')) %>%
  gather('VariableID','val',-MGSampleID) %>%
  mutate(MGSampleID = str_remove_all(MGSampleID, 'M0_')) %>%
  filter(str_detect(MGSampleID, 
                    paste(meta$SampleID, collapse = '|'))) %>%
  left_join(var.names) %>%
  select(-VariableID) %>%
  spread('MGSampleID','val') %>%
  rename(feature = 'DisplayName') %>%
  filter(feature != '-1')

# from cluster, appear to actually be counts but #id_fragment names only
# NOT anonymized !
# motus.hpc <- here('data','metacardis-2021',
#             'metacardis_vs_mOTU_v1_padded_raw_smart_shared_reads_2017-03-26.tsv') %>%
#   read_tsv()
# # BMIs subcohort -- quantitative profiles ??
bmis.feat <- read_tsv(here('data','qmp_genera_abundances_BMIS.tsv'))
bmis.meta <- read_tsv(here('data','metacardis-statins-metadata-raes.tsv')) %>%
  filter(str_detect(SampleID, 'MC')) %>%
  filter(SampleID %in% bmis.feat$SampleID)

table(bmis.meta$`Statin intake`, bmis.meta$Enterotype, useNA='ifany')


## METACARDIS: T2D -- look at drug effects and alcohol ----
discretize <- function(v) {
  q <- quantile(v, probs = seq(0, 1, 0.1), na.rm = TRUE)
  temp <- cut(v, unique(q), include.lowest = TRUE)
  as.double(as.character(factor(temp, 
                                labels = seq_along(levels(temp))))) }

motus.meta <- filter(meta, SampleID %in% colnames(motus.genus))
t2d.meta <- motus.meta %>%
  filter(LABEL %in% c('Healthy','T2D')) %>%
  mutate(LABEL = ifelse(LABEL=='Healthy', -1, 1)) %>%
  mutate_at(vars(contains('_C'), ANTIBIOTIC), 
            ~ ifelse(. == 'Yes', 1, 0)) %>%
  select(SampleID, ALCOHOL_G, everything()) %>%
  mutate(ALCOHOL_G = discretize(ALCOHOL_G)) %>%
  column_to_rownames('SampleID') %>%
  mutate_all(as_factor)

t2d.feat <- select(motus.genus, feature, rownames(t2d.meta)) %>%
  column_to_rownames('feature')
t2d.feat <- t2d.feat[which(rowSums(t2d.feat)!=0),]

## variance explained by T2D and major drug classes from metacardis ----
t2d.varexp <- here('rmd','real_confounder_analysis','proc-data',
                           'metacardis-var-exp-T2D.rds')
if (!file.exists(t2d.varexp)) {
  variances <- list()
  for (i in colnames(t2d.meta)) {
    temp <- t2d.meta[[i]]
    names(temp) <- rownames(t2d.meta)
    if (any(is.na(temp))) { temp <- temp[!is.na(temp)] }
    
    variances[[i]] <- vapply(rownames(t2d.feat), FUN=function(x) {
      x <- t2d.feat[x,names(temp)]
      x <- rank(x)/length(x)
      ss.tot <- sum((x - mean(x))^2)/length(x)
      ss.o.i <- sum(vapply(levels(temp), function(s){
        sum((x[temp==s] - mean(x[temp==s]))^2) }, 
        FUN.VALUE = double(1)))/length(x)
      return(1-ss.o.i/ss.tot) }, FUN.VALUE = double(1)) }
  saveRDS(map_dfr(variances, bind_rows, .id = 'meta.var'),
          here('rmd','real_confounder_analysis','proc-data',
               'metacardis-var-exp-T2D.rds')) 
  variances <- map_dfr(variances, bind_rows, .id = 'meta.var')
} else { 
  variances <- readRDS(here('rmd','real_confounder_analysis','proc-data',
               'metacardis-var-exp-T2D.rds')) }

t2d.varexp <- variances %>%
  gather('genus','var.explained', -meta.var) %>%
  spread('meta.var','var.explained') %>% 
  gather('meta.var','var.explained',-c(genus, LABEL)) %>%
  rename(META = 'var.explained') %>%
  group_by(meta.var) %>%
  mutate(label.rank = rank(LABEL)) %>%
  mutate(meta.rank = rank(META)) %>%
  filter(!str_detect(meta.var, 'STATUS'))

t2d.alcohol <- t2d.varexp %>%
  filter(meta.var == 'ALCOHOL_G') %>%
  ggplot(aes(x = LABEL, y = META)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = 'grey') +
  # facet_wrap(~ meta.var) +
  xlim(0,0.3) + ylim(0,0.3) +
  theme_bw()

# binary metadata for cliff's delta effect size (or gfc) metric
t2d.meta.bin <- t2d.meta %>%
  select(-c('ALCOHOL_G','ANTIBIOTICS_TOTAL'))

fold.changes <- list()
for (i in colnames(t2d.meta.bin)) {
  fold.changes[[i]] <- vapply(rownames(t2d.feat), FUN=function(x) {
    l1 <- levels(get(i, t2d.meta.bin))[1]
    l2 <- levels(get(i, t2d.meta.bin))[2]
    a <- t2d.feat[x,which(get(i, t2d.meta.bin)==l1)]
    b <- t2d.feat[x,which(get(i, t2d.meta.bin)==l2)]
    q.p <- quantile(log10(a + 1e-05), probs = seq(.1, .9, .05), na.rm = TRUE)
    q.n <- quantile(log10(b + 1e-05), probs = seq(.1, .9, .05), na.rm = TRUE)
    return(mean(q.p-q.n)) }, FUN.VALUE = double(1))
}

## cliff's delta for T2D vs major drug classes in metacardis ----
t2d.cliffs <- here('rmd','real_confounder_analysis','proc-data',
                   'metacardis-cliffs-delta-T2D.rds')
if (!file.exists(t2d.cliffs)) {
  cliffs <- list()
  for (i in colnames(t2d.meta.bin)) {
    cliffs[[i]] <- vapply(rownames(t2d.feat), FUN=function(x) {
      l1 <- levels(get(i, t2d.meta.bin))[1]
      l2 <- levels(get(i, t2d.meta.bin))[2]
      a <- t2d.feat[x,which(get(i, t2d.meta.bin)==l1)]
      b <- t2d.feat[x,which(get(i, t2d.meta.bin)==l2)]
      # prob 1>0 or 1>-1 (enriched/depleted in cases/pos classes)
      delta <- orddom::orddom(a, b, onetailed = TRUE)
      return(as.double(delta['delta','ordinal'])) }, 
      FUN.VALUE = double(1)) }
  saveRDS(map_dfr(cliffs, bind_rows, .id = 'meta.var'),
          here('rmd','real_confounder_analysis','proc-data',
               'metacardis-cliffs-delta-T2D.rds'))
  cliffs <- map_dfr(cliffs, bind_rows, .id = 'meta.var')
} else {
  cliffs <- readRDS(here('rmd','real_confounder_analysis','proc-data',
                         'metacardis-cliffs-delta-T2D.rds')) }

# mypal <- ggsci::pal_jama('default', alpha = 0.7)(5)
t2d.cliffs <- cliffs %>%
  gather('genus','delta', -meta.var) %>%
  spread('meta.var','delta') %>% 
  gather('meta.var','delta',-c(genus, LABEL)) %>%
  rename(META = 'delta') %>%
  group_by(meta.var) %>%
  mutate(label.rank = rank(abs(LABEL))) %>%
  mutate(meta.rank = rank(abs(META))) %>%
  filter(!str_detect(meta.var, paste(c('STATUS','RELATED'),
                                     collapse = '|'))) %>%
  mutate(meta.var = str_remove_all(meta.var, '_C')) %>%
  mutate(color = case_when(abs(LABEL) > abs(META) ~ 'ES(label) > ES(drug)', 
                           abs(META) > abs(LABEL) ~ 'ES(drug) > ES(label)', 
                           abs(META) == abs(LABEL) ~ 'ES(drug) == ES(label)',
                           TRUE ~ 'NA'))
ppi.lab <- t2d.cliffs %>% 
  #filter(meta.var == 'PPI') %>% 
  ungroup() %>%
  filter(color == 'ES(drug) > ES(label)') %>% 
  arrange(desc(abs(META)), desc(meta.rank)) %>%
  slice(1:15) 

## cliff's delta for T2D and drugs
t2d.drugs <- t2d.cliffs %>%
  ggplot(aes(x = LABEL, y = META)) +
  geom_abline(slope = 0, intercept = 0, color = 'grey') +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_point(aes(color = color, fill = color)) +
  facet_wrap(~ meta.var) +
  xlim(-0.6, 0.6) + ylim(-0.6, 0.6) +
  labs(x = "Cliff's Delta for T2D", y = "Cliff's Delta for Drug", 
       color = '', fill = '') +
  ggrepel::geom_label_repel(ppi.lab,
                            mapping = aes(label = genus), size = 9/.pt,
                            min.segment.length = unit(0, 'lines'), force = 5) +
  theme_bw() +
  theme(legend.position = 'top', legend.justification = 'left')
ggsave(plot = t2d.drugs, filename = '~/Desktop/metacardis-t2d-drugs.pdf',
       width = 9, height = 7, units = 'in', device = 'pdf')

### METACARDIS: BMIS/OBESITY/METSYN -- PPIs, enterotypes, statins -----
bmi.meta <- bmis.meta %>%
  arrange(SampleID) %>%
  left_join(motus.meta) %>%
  select(1,3,7,16:27, -DIABSTATUS_C, -METFORMIN_C, -DPPIV_C,
         -contains('RELATED')) %>%
  mutate_at(vars(contains('_C'), ANTIBIOTIC), 
            ~ ifelse(. == 'Yes', 1, 0)) %>%
  rename(ET = 'Enterotype') %>%
  rename(BRISTOL = 'Bristol stool score') %>%
  mutate(ALCOHOL_G = discretize(ALCOHOL_G)) %>%
  mutate(ANTIBIOTICS_TOTAL = discretize(ANTIBIOTICS_TOTAL)) %>%
  mutate(LABEL = fct_drop(LABEL)) %>%
  column_to_rownames('SampleID') %>%
  mutate_all(as_factor)

# use the regular genus relative abundances -- not quant. profiles
bmi.feat <- motus.genus %>% 
  select(feature, rownames(bmi.meta)) %>%
  column_to_rownames('feature')
bmi.feat <- bmi.feat[which(rowSums(bmi.feat)!=0),]

## variance explained by obesity/metsyn/drugs/enterotypes ----
bmi.varexp <- here('rmd','real_confounder_analysis','proc-data',
                   'metacardis-var-exp-BMIS.rds')
if (!file.exists(bmi.varexp)) {
  variances <- list()
  for (i in colnames(bmi.meta)) {
    temp <- bmi.meta[[i]]
    names(temp) <- rownames(bmi.meta)
    if (any(is.na(temp))) { temp <- temp[!is.na(temp)] }
    
    variances[[i]] <- vapply(rownames(bmi.feat), FUN=function(x) {
      x <- bmi.feat[x,names(temp)]
      x <- rank(x)/length(x)
      ss.tot <- sum((x - mean(x))^2)/length(x)
      ss.o.i <- sum(vapply(levels(temp), function(s){
        sum((x[temp==s] - mean(x[temp==s]))^2) }, 
        FUN.VALUE = double(1)))/length(x)
      return(1-ss.o.i/ss.tot) }, FUN.VALUE = double(1)) }
  saveRDS(map_dfr(variances, bind_rows, .id = 'meta.var'),
          here('rmd','real_confounder_analysis','proc-data',
               'metacardis-var-exp-BMIS.rds')) 
  variances <- map_dfr(variances, bind_rows, .id = 'meta.var')
} else { 
  variances <- readRDS(here('rmd','real_confounder_analysis','proc-data',
               'metacardis-var-exp-BMIS.rds')) }

bmi.varexp <- variances %>%
  gather('genus','var.explained', -meta.var) %>%
  spread('meta.var','var.explained') %>% 
  gather('meta.var','var.explained',-c(genus, LABEL)) %>%
  rename(META = 'var.explained') %>%
  group_by(meta.var) %>%
  mutate(label.rank = rank(LABEL)) %>%
  mutate(meta.rank = rank(META)) %>%
  mutate(color = case_when(META > LABEL ~ 'var(meta) > var(label)',
                           LABEL > META ~ 'var(label) > var(meta)'))

labeled <- bmi.varexp %>%
  ungroup() %>%
  filter(!str_detect(meta.var, 'ET')) %>%
  arrange(desc(META)) %>%
  slice(1:15)

bmi.metadrugs <- bmi.varexp %>%
  ggplot(aes(x = LABEL, y = META)) +
  geom_abline(slope = 1, intercept = 0, color = 'grey') +
  geom_point(aes(color=color, fill=color)) +
  facet_wrap(~ meta.var) +
  xlim(0,0.3) + ylim(0,0.3) +
  theme_bw() +
  labs(x = "Variance Explained by Obesity/Metabolic Syndrome", 
       y = "Variance Explained by Metavariable", 
       color = '', fill = '') +
  ggrepel::geom_label_repel(labeled,
                            mapping = aes(label = genus), size = 9/.pt,
                            min.segment.length = unit(0, 'lines'), force = 5) +
  theme_bw() +
  theme(legend.position = 'top', legend.justification = 'left')

## cliff's delta for obesity/metsyn/drugs/enterotypes
bmi.meta.bin <- bmis.meta %>%
  arrange(SampleID) %>%
  left_join(motus.meta) %>%
  select(1,3,7,16:27, -ANTIBIOTICS_TOTAL, -DIABSTATUS_C,
         -METFORMIN_C, -DPPIV_C, -contains('RELATED')) %>%
  mutate(L_OBESE = case_when(str_detect(LABEL, 'Healthy') ~ 0,
                             TRUE ~ 1)) %>%
  mutate(L_METSYN = case_when(str_detect(LABEL, 'Healthy') ~ 0,
                              TRUE ~ 1)) %>%
  mutate_at(vars(contains('_C'), ANTIBIOTIC), 
            ~ ifelse(. == 'Yes', 1, 0)) %>%
  rename(ET = 'Enterotype') %>%
  rename(BRISTOL = 'Bristol stool score') %>%
  mutate(ALCOHOL_G = discretize(ALCOHOL_G)) %>%
  mutate(LABEL = fct_drop(LABEL)) %>%
  fastDummies::dummy_cols(c('ET','BRISTOL','ALCOHOL_G','LABEL')) %>%
  select(-c('ET','BRISTOL','LABEL','ALCOHOL_G'), -contains('NA')) %>%
  column_to_rownames('SampleID') %>%
  mutate_all(as.factor)

bmi.cliffs <- here('rmd','real_confounder_analysis','proc-data',
                   'metacardis-cliffs-delta-BMIS.rds')
if (!file.exists(bmi.cliffs)) {
  cliffs <- list()
  for (i in colnames(bmi.meta.bin)) {
    cliffs[[i]] <- vapply(rownames(bmi.feat), FUN=function(x) {
      l1 <- levels(get(i, bmi.meta.bin))[1]
      l2 <- levels(get(i, bmi.meta.bin))[2]
      a <- bmi.feat[x,which(get(i, bmi.meta.bin)==l1)]
      b <- bmi.feat[x,which(get(i, bmi.meta.bin)==l2)]
      # prob 1>0 or 1>-1 (enriched/depleted in cases/pos classes)
      delta <- orddom::orddom(a, b, onetailed = TRUE)
      return(as.double(delta['delta','ordinal'])) }, 
      FUN.VALUE = double(1)) }
  saveRDS(map_dfr(cliffs, bind_rows, .id = 'meta.var'),
          here('rmd','real_confounder_analysis','proc-data',
               'metacardis-cliffs-delta-BMIS.rds'))
  cliffs <- map_dfr(cliffs, bind_rows, .id = 'meta.var')
} else {
  cliffs <- readRDS(here('rmd','real_confounder_analysis','proc-data',
                         'metacardis-cliffs-delta-BMIS.rds')) }

obese.cliffs <- cliffs %>%
  gather('genus','delta', -meta.var) %>%
  spread('meta.var','delta') %>% 
  gather('meta.var','delta',-c(genus, L_OBESE)) %>%
  rename(META = 'delta') %>%
  group_by(meta.var) %>%
  mutate(label.rank = rank(abs(L_OBESE))) %>%
  mutate(meta.rank = rank(abs(META))) %>%
  filter(!str_detect(meta.var, paste(c('STATUS','RELATED'),
                                     collapse = '|'))) %>%
  mutate(meta.var = str_remove_all(meta.var, '_C')) %>%
  mutate(color = case_when(abs(L_OBESE) > abs(META) ~ 'ES(label) > ES(meta)', 
                           abs(META) > abs(L_OBESE) ~ 'ES(meta) > ES(label)', 
                           abs(META) == abs(L_OBESE) ~ 'ES(meta) == ES(label)',
                           TRUE ~ 'NA')) 

labeled <- obese.cliffs %>% 
  ungroup() %>%
  filter(!str_detect(meta.var,'ET_|LABEL')) %>%
  filter(color == 'ES(meta) > ES(label)') %>% 
  arrange(desc(abs(META)), desc(meta.rank)) %>%
  slice(1:15) 

## cliff's delta for T2D and drugs
obese.drugs <- obese.cliffs %>%
  filter(!str_detect(meta.var, 'LABEL|METSYN')) %>%
  filter(!(meta.var %in% paste0('BRISTOL_', seq(2,6,1)))) %>%
  ggplot(aes(x = L_OBESE, y = META)) +
  geom_abline(slope = 0, intercept = 0, color = 'grey') +
  geom_vline(xintercept = 0, color = 'grey') +
  geom_point(aes(color = color, fill = color)) +
  facet_wrap(~ meta.var) +
  xlim(-1, 1) + ylim(-1, 1) +
  labs(x = "Cliff's Delta for Obesity", y = "Cliff's Delta for Metavariable", 
       color = '', fill = '') +
  ggrepel::geom_label_repel(labeled,
                            mapping = aes(label = genus), size = 9/.pt,
                            min.segment.length = unit(0, 'lines'), force = 5) +
  theme_bw() +
  theme(legend.position = 'top', legend.justification = 'left')
ggsave(plot = obese.drugs, filename = '~/Desktop/metacardis-obesity-drugs.pdf',
       width = 9, height = 7, units = 'in', device = 'pdf')