# ##############################################################################
#
## Create simulations with the negative binominal method from Hawinkel et al.
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
simulation <- paste0(temp.loc, 'sim_negbin_cor.h5')
sim.method <- 'negbin'
rep <- 10

# create simulations
simulate.data(feat = feat.zeevi, meta=meta.zeevi,
              sim.location = simulation,
              sim.type='cross-section',
              sim.method=sim.method,
              filt.params = list(ab.cutoff=as.numeric(params$ab.cutoff),
                                 prev.cutoff=params$prev.cutoff,
                                 log.n0=as.numeric(params$log.n0)),
              sim.params = list(ab.scale=params$ab.scale,
                                prop.markers=params$prop.markers,
                                class.balance=params$class.balance,
                                correlation=TRUE,
                                repeats=rep))

# create testing indices
create.test.idx(sim.location=simulation,
                subsets=params$subset_size,
                repetitions=50)

# move file to permament place
file.move(simulation, here('simulations', 'others'))

