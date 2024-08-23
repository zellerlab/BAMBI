# ##############################################################################
#
## Create simulations with implantation method
##      feature implantation with Study as confounder
##      feature implantation with an artificial confounder
##      feature implantation with a global confounder (?)
#
# ##############################################################################

library("here")
library("tidyverse")
library("filesstrings")

# load the data
source(here("create_simulations", "load_data.R"))
dataset.list$Zeevi_KEGG <- NULL
devtools::load_all(simba.loc)

params <- yaml::yaml.load_file(here('create_simulations', "parameters.yaml"))
# temp.loc definition missing?

# parameters
sim.method <- 'resampling'
rep <- 50
conf.bias <- 0.5

# ##############################################################################
# Study confounder
# simulate and add both sets of test indices

combinations <- c('Zeevi_WGS-Schirmer_WGS', 'Zeevi_WGS-TwinsUK_WGS',
                  'Zeevi_WGS-MI_WGS', 'Schirmer_WGS-TwinsUK_WGS',
                  'Schirmer_WGS-MI_WGS', 'TwinsUK_WGS-MI_WGS')
for (i in combinations){
  studies <- str_split(i, pattern = '-')[[1]]
  simulation <- paste0(temp.loc, 'sim_study_conf_', studies[1], '-', 
                       studies[2], '.h5')
  message(simulation)
  if (file.exists(here('simulations', 'conf_sim', 
                       paste0('sim_study_conf_', studies[1], '-', 
                              studies[2], '.h5')))){next()}
  # create simulations
  create.data.simulation(feat = list(dataset.list[[studies[1]]]$feat, 
                                     dataset.list[[studies[2]]]$feat),
                meta = list(dataset.list[[studies[1]]]$meta, 
                            dataset.list[[studies[2]]]$meta),
                sim.location = simulation,
                sim.type='cross-section',
                sim.method=sim.method,
                filt.params = list(ab.cutoff=as.numeric(params$ab.cutoff),
                                   prev.cutoff=params$prev.cutoff,
                                   log.n0=as.numeric(params$log.n0)),
                sim.params = list(ab.scale=params$ab.scale,
                                  prev.scale=params$prev.scale,
                                  prop.markers=params$prop.markers,
                                  class.balance=params$class.balance,
                                  conf='batch',
                                  conf.params=list(bias=conf.bias, prop=0.5),
                                  prop.markers=params$prop.markers,
                                  feature.type='all',
                                  repeats=rep))
  message("Created study-conf simulation!")

  # create normal testing indices
  create.test.idx(sim.location=simulation,
                  subsets=c(50, 100, 200),
                  repetitions=50)
  message("Created norma test indices for study-conf simulation!")

  # create resampled test indices
  create.biased.idx(sim.location=simulation,
                    lab.svar = NULL, meta.svar = NULL, 
                    subsets = c(50, 100, 200, 400),
                    repetitions = 50,
                    strat.scheme = 'confounder',
                    bias=c(0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95))
  message("Created confounded test indices for study-conf simulation!")

  # move to permanent place
  file.move(simulation, here('simulations', 'conf_sim'))
}



# ##############################################################################
# Global confounder
# simulate and add both sets of test indices

# simulation <- paste0(temp.loc, 'sim_global_conf')

# create simulations
# create.data.simulation(feat = feat.zeevi,
#               meta = meta.zeevi,
#               sim.location = paste0(simulation, '.h5'),
#               sim.type='cross-section',
#               sim.method=sim.method,
#               filt.params = list(ab.cutoff=as.numeric(params$ab.cutoff),
#                                  prev.cutoff=params$prev.cutoff,
#                                  log.n0=as.numeric(params$log.n0)),
#               sim.params = list(ab.scale=params$ab.scale,
#                                 prev.scale=params$prev.scale,
#                                 prop.markers=params$prop.markers,
#                                 class.balance=params$class.balance,
#                                 conf='global',
#                                 conf.params=list(bias=conf.bias, prop=0.5),
#                                 prop.markers=params$prop.markers,
#                                 feature.type='all',
#                                 repeats=rep))

# create normal testing indices
# create.test.idx(sim.location=paste0(simulation, '.h5'),
#                 subsets=c(50, 100, 200),
#                 repetitions=50)

# create resampled test indices
# create.biased.idx(sim.location=paste0(simulation, '.h5'),
#                   lab.svar = NULL, meta.svar = NULL, 
#                   subsets = c(50, 100, 200, 400),
#                   repetitions = 50,
#                   strat.scheme = 'confounder',
#                   bias=c(0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95))

# move to permanent place
# file.move(paste0(sim, '.h5'), here('simulations', 'conf_sim'))

# ##############################################################################
# Artificial confounder
# simulate and add both sets of test indices
conf.feat.prop <- c(0.05, 0.1, 0.2, 0.5)
datasets <- names(dataset.list)[str_detect(names(dataset.list), 'WGS')]
for (cfp in conf.feat.prop){
  for (i in names(dataset.list)){
      simulation <- paste0(temp.loc, 'sim_artificial_conf_', cfp, '_', i, '.h5')
      message(simulation)
      if (file.exists(here('simulations', 'conf_sim', 
                           paste0('sim_artificial_conf_', cfp, 
                                  '_', i, '.h5')))){next()}
    
      # create simulations
      create.data.simulation(feat = dataset.list[[i]]$feat,
                    meta = dataset.list[[i]]$meta,
                    sim.location = simulation,
                    sim.type='cross-section',
                    sim.method=sim.method,
                    filt.params = list(ab.cutoff=as.numeric(params$ab.cutoff),
                                       prev.cutoff=params$prev.cutoff,
                                       log.n0=as.numeric(params$log.n0)),
                    sim.params = list(ab.scale=params$ab.scale,
                                      prev.scale=params$prev.scale,
                                      prop.markers=params$prop.markers,
                                      class.balance=params$class.balance,
                                      conf='artificial',
                                      conf.params=list(bias=conf.bias, 
                                                       prop=0.5,
                                                       feat.prop=cfp),
                                      prop.markers=params$prop.markers,
                                      feature.type='all',
                                      repeats=rep))
    
      # create normal testing indices
      create.test.idx(sim.location=simulation,
                      subsets=c(50, 100, 200),
                      repetitions=50)
    
      # create resampled test indices
      create.biased.idx(sim.location=simulation,
                        lab.svar = NULL, meta.svar = NULL, 
                        subsets = c(50, 100, 200, 400),
                        repetitions = 50,
                        strat.scheme = 'confounder',
                        bias=c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9))
      
      # move to permanent place
      file.move(simulation, here('simulations', 'conf_sim'))
  }
}
