# ##############################################################################
#
## Scripts to check the reality of simulations using batchtools
#
# ##############################################################################

library("tidyverse")
library("batchtools")
library("devtools")
library("SIAMCAT")
library("here")

args <- commandArgs(trailingOnly = TRUE)

# load package
devtools::load_all(simba.loc)

# get parameters
simulation <- args[1]

sim.files <- list.files(here('simulations'), 
                        full.names = TRUE, recursive = TRUE)
sim.file <- sim.files[str_detect(sim.files, pattern = simulation)]
if (length(sim.file) != 1){
  stop("Cannot find the correct simulation!")
}
if (!file.exists(sim.file)){
  stop("No such simulation exists!")
}

# copy simulation file to scratch for faster I/O access
file.loc <- paste0(temp.loc, 'simulations/', simulation, '.h5')
file.copy(from=sim.file, file.loc)

# start of real script :D
sim.params <- h5read(file.loc, name='simulation_parameters/sim_params')

# get groups
temp <- h5ls(file.loc, recursive=FALSE)
groups <- temp$name
groups <- setdiff(groups, c('original_data', 'simulation_parameters'))

groups <- unique(str_remove(groups, '_rep[0-9]*$'))

# create a job registry
job.registry <- here(temp.loc, 'reality_checks_registries', simulation)
if (!dir.exists(job.registry)){
    makeRegistry(file.dir=job.registry, work.dir=here(),
    conf.file=here('cluster_config', 'batchtools_reality_conf.R'))
} else {
    loadRegistry(file.dir=job.registry, writeable=TRUE,
                 conf.file=here('cluster_config', 'batchtools_reality_conf.R'))
}


.f <- function(g, sim.file, max.iter, simba.loc){
    # load SIMBA
    devtools::load_all(simba.loc)

    res.list <- list()

    # iterate over repetitions
    for (x in seq_len(max.iter)){
        res.list[[x]] <- reality.check(sim.location=sim.file,
                                       group=paste0(g, '_rep', x),
                                       ml=TRUE)
    }
    return(res.list)
}

# add jobs
batchMap(.f, groups, more.args=list('sim.file'=sim.file,
                                    'max.iter'=sim.params$repeats, 
                                    'simba.loc'=simba.loc))
submitJobs()
