# ##############################################################################
#
## Evaluate all test results for a single simulation
#
# ##############################################################################

library("tidyverse")
library("batchtools")
library("here")

# load SIMBA
devtools::load_all(simba.loc)

args <- commandArgs(trailingOnly = TRUE)
simulation <- args[1]
test <- args[2]
if (!is.na(test)){
  message("Evaluating jobs for the test: ", test)
} else {
  message("Evaluating all different tests!")
}

# get parameters
sim.files <- list.files(here('simulations'), 
    full.names = TRUE, recursive = TRUE)
sim.file <- sim.files[str_detect(sim.files, pattern = simulation)]
if (length(sim.file) != 1){
    stop("Cannot find the correct simulation!")
}
if (!file.exists(sim.file)){
    message('Evaluating single simulation!')
}

# load source registry for simulation
job.registry <- here('test_results_registries', paste0(simulation, '_test'))
if (!dir.exists(job.registry)){
 stop("Simulation hasn't been tested yet! Simulations which have include ",
       paste(setdiff(list.files(here('test_results_registries')),
                     c('README.md')), collapse = ', '))
} else {
  test.reg <- loadRegistry(file.dir = job.registry,
                           conf.file = here('cluster_config',
                                            'batchtools_test_conf.R'),
                           make.default = TRUE)
}

# get job parameters
parameters <- getJobTable() %>%
  as_tibble() %>%
  select(job.id, done, error, time.running, mem.used, problem,
         algorithm, algo.pars) %>%
  unnest_wider(algo.pars)
if (!is.na(test)){
  parameters <- parameters %>%
    filter(algorithm==!!test)
}

# find completed jobs
parameters <- parameters %>%
  mutate(done.bool = !is.na(done) & is.na(error))

# print a statement about the status?
if (!all(parameters$done.bool)){
  print("Not all jobs have finished running yet!")
  print("Missing jobs are for the following test: ")
  print(parameters %>% filter(!done.bool) %>% group_by(test) %>% tally())
}

# result file location
if (!is.na(test)){
  fn.res <- here('test_results', paste0(simulation, '_', test, '.tsv'))
} else {
  fn.res <- here('test_results', paste0(simulation, '.tsv'))
}

# check if some things have been evaluated before
if (file.exists(fn.res)){
  temp.res <- read_tsv(fn.res)
  parameters <- parameters %>% 
    filter(!job.id %in% temp.res$job.id)
  message('Initial results already exists! Picking up where we left...')
}

# extract result tables
message("Reducing results ...")
pb <- progress_bar$new(total=nrow(parameters))
# loop through the jobs which are done
for (i in seq_len(nrow(parameters))){
  res.list <- list()
  job.id <- parameters$job.id[i]
  if (!parameters$done.bool[i]){
    next()
  }
  # retrieve results
  res <- loadResult(job.id)
  for (g in names(res)){
    for (s in names(res[[g]])){
      eval.res <- eval.test(sim.file, g, res[[g]][[s]], alpha=0.05)
      res.list[[(length(res.list)+1)]] <- eval.res %>% 
        as_tibble() %>% mutate(group=g, subset=s) %>% 
        mutate(job.id=job.id)
    }
  }
  res.list <- bind_rows(res.list)
  
  # bring into right shape
  tmp <- res.list %>% 
    mutate(problem=parameters$problem[i],
           test=parameters$test[i],
           norm=parameters$norm[i],
           time.running=parameters$time.running[i],
           mem.used=parameters$mem.used[i])
  # save
  if (!file.exists(fn.res)){
    write_tsv(tmp, file = fn.res)
  } else {
    write_tsv(tmp, file = fn.res, append = TRUE)
  }
  pb$tick()
}

message("Evaluation results saved in test_results!")
