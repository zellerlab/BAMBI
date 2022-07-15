# ##############################################################################
#
## Create simulations with the beta binominal method from Hawinkel et al.
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
simulation <- paste0(temp.loc, 'sim_betabin.h5')
sim.method <- 'betabin'
rep <- 10

message(simulation)
message(sim.method)

# create simulations
create.data.simulation(feat = feat.zeevi, meta=meta.zeevi,
              sim.location = simulation,
              sim.type='cross-section',
              sim.method=sim.method,
              filt.params = list(ab.cutoff=as.numeric(params$ab.cutoff),
                                 prev.cutoff=params$prev.cutoff,
                                 log.n0=as.numeric(params$log.n0)),
              sim.params = list(ab.scale=params$ab.scale,
                                prop.markers=params$prop.markers,
                                class.balance=params$class.balance,
                                repeats=rep))

# create testing indices
create.test.idx(sim.location=simulation,
                subsets=params$subset_size,
                repetitions=50)

# move file to permament place
file.move(simulation, here('simulations', 'others'))

