# extra metadata columns needed for SIMBA
genus.map <- here('data','metacardis',
                   'hub.taxon.motus.Genus.raw.mgsamples.v3.r') %>%
  read_tsv() %>%
  select(MGSampleID) %>%
  separate(MGSampleID, '_', into = c('Timepoint','Individual_ID')) %>%
  mutate(Timepoint = case_when(Timepoint == 'M0' ~ 1,
                               Timepoint == 'M3' ~ 2,
                               Timepoint == 'M6' ~ 3,
                               Timepoint == 'M12' ~ 4))

# not needed, all present in genus data
# species.map <- here('data','metacardis',
#                     'hub.cellcount.motu.Species.mgsamples.v2.r') %>%
#   read_tsv() %>%
#   select(MGSampleID) %>%
#   separate(MGSampleID, '_', into = c('Timepoint','Individual_ID')) %>%
#   mutate(Timepoint = case_when(Timepoint == 'M0' ~ 1,
#                                Timepoint == 'M3' ~ 2,
#                                Timepoint == 'M6' ~ 3,
#                                Timepoint == 'M12' ~ 4))
# sample.map <- distinct(bind_rows(species.map, genus.map))

meta.vars <- here('data','metacardis', 'hub.metadata.variables.v10.r') %>%
  read_tsv() %>%
  filter(str_detect(Category, 'Drugs')) %>% 
  filter(VariableType == 'Categorical')
motus.meta <- here('data','metacardis', 'hub.metadata.samples.v10.r') %>%
  read_tsv() %>%
  mutate(Label = fct_collapse(PATGROUPFINAL_C,
                              # healthy controls
                              HC = c('8','1 RH'),
                              # metabolic syndrome, no CCAD/T2D
                              MS = '1',
                              # severe obesity (>35 bmi), no CCAD/T2D
                              OB = '2a',
                              # candidates for bariatric surgery (comorbid, >40bmi)
                              OB_BAR = '2b',
                              # T2D no CCAD
                              T2D = '3',
                              # first coronary artery disease event
                              CAD = '4',
                              # chronic coronary artery disease, no heart failure
                              CCAD = '5',
                              # chronic CAD with heart failure
                              CCAD_HF = '6',
                              # heart failure, not from CAD
                              HF = '7')) %>%
  select(SampleID, Label, all_of(meta.vars$VariableID)) %>%
  mutate(across(c('ANTIBIOTIC', contains('_C')), 
                ~ fct_recode(., `1` = 'Yes', `0` = 'No'))) %>%
  rename_with(tolower, contains('_C')) %>%
  rename_with(~ str_remove_all(., '_c')) %>%
  rename(antibiotics = 'ANTIBIOTIC') %>%
  left_join(genus.map, by = c('SampleID' = 'Individual_ID')) %>%
  rename(Individual_ID = 'SampleID') %>%
  # mutate(Sample_ID = paste0(Individual_ID, '_', Timepoint)) %>%
  mutate(Sample_ID = Individual_ID) %>%
  select(contains('_ID'), Timepoint, Label, everything()) %>%
  filter(Timepoint == 1)

# raw, genus-level relative abundances
gen.vars <- here('data','metacardis',
                  'hub.taxon.motus.Genus.raw.variables.v3.r') %>%
  read_tsv() %>% select(VariableID, DisplayName)
motus.genus <- here('data','metacardis',
                    'hub.taxon.motus.Genus.raw.mgsamples.v3.r') %>%
  read_tsv() %>%
  separate(MGSampleID, '_', into = c('Timepoint','Individual_ID')) %>%
  filter(Timepoint == 'M0') %>%
  gather('VariableID','val', -c(Individual_ID, Timepoint)) %>%
  filter(Individual_ID %in% motus.meta$Individual_ID) %>%
  left_join(gen.vars) %>%
  select(-VariableID, -Timepoint) %>%
  # filter(DisplayName != '-1') %>%
  spread('Individual_ID','val') %>%
  rename(feature = 'DisplayName') %>%
  filter(rowSums(.[,2:ncol(.)]) > 0) %>%
  column_to_rownames('feature')

# cell count corrected abundances (used for 2019 simulations)
# spec.vars <- here('data','metacardis',
#                  'hub.cellcount.motu.Species.variables.v2.r') %>%
#   read_tsv() %>%
#   rename(VariableID = 'Variable ID', DisplayName = 'Display name') %>% 
#   select(VariableID, DisplayName)
# contains c. bartlettii
motus.species <- here('data','metacardis',
                      'hub.cellcount.motu.Species.mgsamples.v2.r') %>%
  read_tsv() %>%
  separate(MGSampleID, '_', into = c('Timepoint','Individual_ID')) %>%
  filter(Timepoint == 'M0') %>%
  gather('feature','val', -c(Individual_ID, Timepoint)) %>%
  filter(Individual_ID %in% motus.meta$Individual_ID) %>%
  select(-Timepoint) %>%
  # filter(feature != '-1') %>%
  filter(!str_detect(feature, 'linkage')) %>%
  mutate(feature = str_replace_all(feature, '__', '_')) %>% 
  mutate(feature = str_remove_all(feature, '^[_]')) %>%
  spread('Individual_ID','val') %>%
  filter(rowSums(.[,2:ncol(.)]) > 0) %>%
  column_to_rownames('feature')

rm(gen.vars, meta.vars, genus.map)