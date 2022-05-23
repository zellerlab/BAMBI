# ##############################################################################
#
## Resample and permute real data and save in .h5 file like simulated data 
##      metacardis
#
# ##############################################################################

library('here')
library('tidyverse')
simba.loc <- '~/Dropbox/_phd/projects/sim-bench/collab-jakob/SIMBA'
devtools::load_all(simba.loc)

# load the data -- don't combine with zeevi/twins
source(here('create_simulations', 'load_metacardis.R'))

params <- yaml::yaml.load_file(here('create_simulations', 'parameters.yaml'))

diseases <- c('T2D','CCAD','OB','MS')
tax.res <- 'species'

for (phenotype in diseases) {
  
  meta <- motus.meta %>%
    filter(str_detect(Label, phenotype) | Label == 'HC') %>%
    mutate(Label = fct_drop(Label)) %>%
    fastDummies::dummy_cols(select_columns = 'Label',
                            remove_selected_columns = TRUE) %>%
    select(contains(c('_ID','Time','Label')), everything()) %>%
    mutate_if(is.factor, ~ as.numeric(as.character(.))) %>%
    column_to_rownames('Sample_ID')
  if (tax.res == 'species') {
    meta <- filter(meta, rownames(meta) %in% colnames(motus.species))
    df.feat <- motus.species[,rownames(meta)]
  } else if (tax.res == 'genus') { 
    df.feat <- motus.genus[,rownames(meta)] 
  } else { stop('taxonomic level not recognized!') }
  
  # dummy cols not needed for T2D
  if (phenotype == 'T2D' | phenotype == 'MS') 
    meta <- select(meta, -Label_HC) 
  
  # possible ways to pool phenotypes
  labels <- meta %>%
    select(contains('Label')) %>%
    colnames() %>%
    str_remove_all('Label_')
  
  for (config in labels) {
    # can never be 1/1 overlap with dummy vars; remove extra labels
    df.meta <- meta %>%
      rename(tmp = paste0('Label_', config)) %>%
      select(-contains('Label')) %>%
      rename(Label = 'tmp') 
    # default -- pool disease signal vs healthy controls
    if (config == 'HC') {
      df.meta <- df.meta %>%
        mutate(Label = ifelse(Label == 1, -1, 1))
      simulation <- here('simulations','conf_real',
                         paste0('metacardis_', tax.res, '_drugs_', phenotype,
                                '_pooled'))
    } else {
      df.meta <- df.meta %>%
        mutate(Label = ifelse(Label == 1, 1, -1))
      simulation <- here('simulations','conf_real',
                         paste0('metacardis_', tax.res, '_drugs_', config)) }
    message('+ generating ', simulation)
    
    # create "simulations"
    # sim.type should be data.type
    # sim.method should be sim.type (pass, resampling, parametric) +
    #                      sim.method (NULL, conf/vanilla, e.g. weiss)
    simulate.data(feat = df.feat,
                  meta = df.meta,
                  sim.location = paste0(simulation, '.h5'),
                  sim.type = 'cross-section',
                  sim.method = 'pass',
                  filt.params = list(ab.cutoff=as.numeric(params$ab.cutoff),
                                     prev.cutoff=params$prev.cutoff,
                                     log.n0=as.numeric(params$log.n0)),
                  sim.params = NULL)
    
    # save additional 'real' data and metadata
    if (file.exists(paste0(simulation, '.h5'))) {
      sim <- h5read(paste0(simulation, '.h5'), 'original_data')
      
      # save real biological names
      feat.names <- tibble(original = rownames(df.feat),
                           renamed = sim$feature_names)
      filt.feat.names <- feat.names %>%
        filter(renamed %in% sim$filt_feature_names)
      h5write(obj = feat.names$original,
              file = paste0(simulation, '.h5'),
              name = 'original_data/feature_names_real')
      h5write(obj = filt.feat.names$original,
              file = paste0(simulation, '.h5'),
              name = 'original_data/filt_feature_names_real')
      h5write(obj = rownames(df.meta),
              file = paste0(simulation, '.h5'),
              name = 'original_data/sample_names_real') 
      } else { stop('Simulation not created yet!') }
    
    # get naive idx for each subset size, saves @ /original_data 
    create.test.idx(sim.location = paste0(simulation, '.h5'),
                    subsets = c(50, 100, 200),
                    repetitions = 50)
    
    create.biased.idx(sim.location = paste0(simulation, '.h5'),
                      lab.svar = 'Label',
                      meta.svar = setdiff(colnames(df.meta),
                                          c('Individual_ID','Timepoint')),
                      subsets = c(50, 100, 200),
                      repetitions = 50,
                      strat.scheme = 'random')
    
    create.biased.idx(sim.location = paste0(simulation, '.h5'),
                      lab.svar = 'Label',
                      meta.svar = setdiff(colnames(df.meta),
                                          c('Individual_ID','Timepoint')),
                      subsets = c(50, 100, 200),
                      repetitions = 50,
                      strat.scheme = 'confounder',
                      bias = c(0.5, 0.75))
  }
}