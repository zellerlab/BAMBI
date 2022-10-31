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

# extend for the different adjustments
parameters <- parameters %>% 
  select(job.id, done.bool, test, norm, time.running, mem.used, problem) %>% 
  mutate(adjust1='none', adjust2='BH', adjust3='BY') %>% 
  pivot_longer(cols=c(adjust1, adjust2, adjust3),
               names_to = 'name', values_to = 'adjust') %>% 
  select(-name) %>% 
  mutate(eval=paste0(job.id, '-', adjust)) %>% 
  filter(done.bool)

# check if some things have been evaluated before
if (file.exists(fn.res)){
  temp.res <- read_tsv(fn.res) %>% 
    mutate(eval=paste0(job.id, '-', adjust))
  parameters <- parameters %>% 
    filter(!eval %in% temp.res$eval)
  message('Initial results already exists! Picking up where we left...')
  message('Number of missing evals: ', nrow(parameters))
}

# extract result tables
message("Reducing results ...")
pb <- progress_bar$new(total=length(unique((parameters$job.id))))
# loop through the jobs which are done
for (i in seq_along(unique(parameters$job.id))){
  res.list <- list()
  # retrieve results
  res <- loadResult(i)
  tmp.param <- parameters %>% 
    filter(job.id==i)
  for (g in names(res)){
    for (s in names(res[[g]])){
      for (x in tmp.param$adjust){
        eval.res <- eval.test(sim.file, g, res[[g]][[s]], adjust=x, alpha=0.05)
        res.list[[(length(res.list)+1)]] <- eval.res %>% 
          as_tibble() %>% mutate(group=g, subset=s, adjust=x) %>% 
          mutate(job.id=tmp.param$job.id)
      }
    }
  }
  res.list <- bind_rows(res.list)
  
  # bring into right shape
  tmp <- res.list %>% 
    mutate(problem=tmp.param$problem[1],
           test=tmp.param$test[1],
           norm=tmp.param$norm[1],
           time.running=tmp.param$time.running[1],
           mem.used=tmp.param$mem.used[1])
  # save
  if (!file.exists(fn.res)){
    write_tsv(tmp, file = fn.res)
  } else {
    write_tsv(tmp, file = fn.res, append = TRUE)
  }
  pb$tick()
}

message("Evaluation results saved in test_results!")
