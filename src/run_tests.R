# ##############################################################################
#
## Run all tests for a single simulation
#
# ##############################################################################


library("tidyverse")
library("batchtools")
library("rhdf5")
library("here")
library("yaml")

# load SIMBA
devtools::load_all(simba.loc)

# get sim & paramaters
args <- commandArgs(trailingOnly = TRUE)
simulation <- args[1]
test <- args[2]

message(simulation)
message(test)

if (!is.na(test)){
  message("Submitting jobs for the test: ", test)
} else {
  message("Submitting all different tests!")
}


# make experiment registry - might need to load pkgs here
# load source registry for simulation
job.registry <- paste0(temp.loc, 'test_results_registries/', 
                       paste0(simulation, '_test'))
message(job.registry)
if (!dir.exists(job.registry)){
  stop("Simulation hasn't been tested yet! Simulations which have include ",
       paste(setdiff(list.files(paste0(temp.loc, 'test_results_registries')),
                     c('README.md')), collapse = ', '))
} else {
  test.reg <- loadRegistry(file.dir = job.registry,
                           conf.file = here('cluster_config',
                                            'batchtools_test_conf.R'),
                           writeable = TRUE)
}

print(getStatus())
errs <- findErrors()
if (nrow(errs) > 0){
	resetJobs(errs)
	print(getStatus())
}

if (is.na(test)){
  not.submitted <- findNotDone()
  print(nrow(not.submitted))
  submitJobs(not.submitted$job.id, resources=list(max.concurrent.jobs=2000))
} else {
  not.submitted <- getJobPars(findNotSubmitted()) %>% filter(algorithm == test)
  print(nrow(not.submitted))
  submitJobs(not.submitted$job.id, resources=list(max.concurrent.jobs=2000))  
}
