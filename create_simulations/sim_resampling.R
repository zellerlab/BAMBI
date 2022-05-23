# ##############################################################################
#
## Create simulations with resampling method
##      feature implantation in different abundance schemes
#
# ##############################################################################

library("here")
library("tidyverse")
library("yaml")
library("devtools")
library("filesstrings")

# load the data
source(here("create_simulations", "load_data.R"))

# SIMBA
devtools::load_all(simba.loc)

# parameters
params <- yaml.load_file(here('create_simulations', "parameters.yaml"))

# parameters
feature.type <- c('all', 'low', 'middle', 'high', 'low_middle', 
                  'abundance', 'inverse-abundance')
simulation <- paste0(temp.loc, 'simulations/sim_')

for (ft in feature.type){
  sim <- paste0(simulation, ft, '.h5')
  
  # create simulations
  simulate.data(feat = feat.zeevi,
                meta = meta.zeevi,
                sim.location = sim,
                sim.type='cross-section',
                sim.method='resampling',
                filt.params = list(ab.cutoff=as.numeric(params$ab.cutoff),
                                   prev.cutoff=params$prev.cutoff,
                                   log.n0=as.numeric(params$log.n0)),
                sim.params = list(ab.scale=params$ab.scale,
                                  prev.scale=params$prev.scale,
                                  prop.markers=params$prop.markers,
                                  conf=NULL,
                                  class.balance=params$class.balance,
                                  feature.type=ft,
                                  repeats=100))
  
  # create testing indices
  create.test.idx(sim.location=sim,
                  subsets=params$subset_size,
                  repetitions=50)
  file.move(sim, here('simulations', 'resampling'))
}


# additionally create un-balanced implantation for compositional effects
sim <- paste0(simulation, 'compositional.h5')

# create simulations
simulate.data(feat = feat.zeevi,
              meta = meta.zeevi,
              sim.location = sim,
              sim.type='cross-section',
              sim.method='resampling',
              filt.params = list(ab.cutoff=as.numeric(params$ab.cutoff),
                                 prev.cutoff=params$prev.cutoff,
                                 log.n0=as.numeric(params$log.n0)),
              sim.params = list(ab.scale=params$ab.scale,
                                prev.scale=params$prev.scale,
                                prop.markers=params$prop.markers,
                                conf=NULL,
                                class.balance=params$class.balance,
                                feature.type='all',
                                balanced=FALSE,
                                repeats=100))
create.test.idx(sim.location=sim,
                subsets=params$subset_size,
                repetitions=50)
file.move(sim, here('simulations', 'resampling'))
