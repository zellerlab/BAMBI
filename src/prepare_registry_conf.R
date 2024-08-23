# ##############################################################################
#
## Run tests for a confounded simulation
#
# ##############################################################################

library("tidyverse")
library("batchtools")
library("here")

# load SIMBA
devtools::load_all(simba.loc)

# get sim & paramaters
args <- commandArgs(trailingOnly = TRUE)
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
file.loc <- paste0(temp.loc, 'simulations/', simulation, '.h5')
# file.copy(from=sim.file, file.loc)

sim.params <- h5read(file.loc, name='simulation_parameters/sim_params')

# tests <- list or vector
tests <- c('limma', 'wilcoxon', 'lm', 'ANCOMBC', 'metagenomeSeq2',
           'wilcoxon_conf', 'limma_conf', 'ANCOMBC_conf',
           'lme_conf', 'fastANCOM', 'fastANCOM_conf', "corncob", "corncob_conf")

# make experiment registry - might need to load pkgs here
job.registry <- paste0(temp.loc, 'test_results_registries/', 
                       paste0(simulation, '_test'))
if (!dir.exists(job.registry)) {
  reg <- makeExperimentRegistry(file.dir = job.registry, 
                                work.dir = here(),
                                conf.file = here('cluster_config',
                                                'batchtools_test_conf.R'),
                                make.default = TRUE) 
} else {
  stop("Registry already exists!")
}

# get groups -- e.g. ab1_rep1
temp <- h5ls(file.loc, recursive=FALSE)
groups <- setdiff(temp$name, c('original_data', 'simulation_parameters'))
conf.params <- h5read(file.loc, name='simulation_parameters/conf_bias_test_idx')
bias.terms <- conf.params$bias
subsets <- conf.params$subsets

### ADD PROBLEMS, ALGORITHMS, THEN EXPERIMENTS
sim.location <- file.loc

# add problems
addProblem(name = simulation,
           data = sim.location,
           fun = function(data, job) data,
           reg = reg)

# add algorithms (tests)
.f_apply <- function(data, job, instance, group, subsets, norm, test, 
                     bias, ...){
  
  type <- paste0('confounder::', bias)
  if (grepl('_conf', test)){
    conf <- 'conf'
    test <- gsub(x=test, pattern='_conf', replacement='')
  } else {
    conf <- NULL
  }
  apply.test(sim.location = data,
             group=group,
             subset=subsets,
             norm=norm,
             test=test,
             type=type,
             conf=conf)
}
walk(tests, ~ addAlgorithm(name = .,  fun = .f_apply))

# a bit manual tweaking, but could work...

ades.fast <- c('fastANCOM', 'fastANCOM_conf') %>%
  map(~ data.table::CJ(test = .,
                       group=unique(str_remove(groups, '_rep[0-9]*$')),
                       norm = c('pass'),
                       subsets=NA_real_,
                       bias=seq_along(bias.terms))) %>%
  set_names(c('fastANCOM', 'fastANCOM_conf'))
addExperiments(algo.designs = ades.fast)
submitJobs()

# corncob
groups <- groups[str_detect(groups, 'rep[1-9]{1}$')]
ades.ancombc <- c('corncob', 'corncob_conf') %>%
  map(~ data.table::CJ(test = .,
                       group=groups,
                       norm = c('pass'),
                       subsets=NA_real_,
                       bias=seq_along(bias.terms))) %>%
  set_names(c('corncob', 'corncob_conf'))

# slow tests
ades.ancombc <- c('ANCOMBC', 'ANCOMBC_conf') %>%
  map(~ data.table::CJ(test = .,
                     group=groups,
                     norm = c('pass'),
                     subsets=NA_real_,
                     bias=seq_along(bias.terms))) %>%
  set_names(c('ANCOMBC', 'ANCOMBC_conf'))

# limma
ades.limma <- c('limma', 'limma_conf') %>%
  map(~ data.table::CJ(test=.,
                       group=unique(str_remove(groups, '_rep[0-9]*$')),
                       norm = c('TSS.arcsin'),
                       subsets=NA_real_,
                       bias=seq_along(bias.terms))) %>%
  set_names(c('limma', 'limma_conf'))

# LM
ades.lm <- c('lm') %>% 
  map(~ data.table::CJ(test=.,
                       group=unique(str_remove(groups, '_rep[0-9]*$')),
                       norm = c('TSS.log'),
                       subsets=NA_real_,
                       bias=seq_along(bias.terms))) %>%
  set_names(c('lm'))

ades.lme.c <- c('lme_conf') %>%
  map(~ data.table::CJ(test=.,
                       group=groups,
                       norm = c('TSS.log'),
                       subsets=NA_real_,
                       bias=seq_along(bias.terms))) %>%
  set_names(c('lme_conf'))

# wilcoxon
ades.wi <- c('wilcoxon') %>%
  map(~ data.table::CJ(test=.,
                       group=unique(str_remove(groups, '_rep[0-9]*$')),
                       norm = c('TSS'),
                       subsets=NA_real_,
                       bias=seq_along(bias.terms))) %>%
  set_names(c('wilcoxon'))

# wilcoxon_conf
ades.wi.c <- c('wilcoxon_conf') %>%
  map(~ data.table::CJ(test=.,
                       group=groups,
                       norm = c('TSS'),
                       subsets=NA_real_,
                       bias=seq_along(bias.terms))) %>%
  set_names(c('wilcoxon_conf'))


# others
ades.m <- c('metagenomeSeq2') %>%
  map(~ data.table::CJ(test = .,
                       group=groups,
                       norm = c('pass'),
                       subsets=NA_real_),
                       bias=seq_along(bias.terms)) %>%
  set_names(c('metagenomeSeq2'))

ades <- append(ades.ancombc, ades.limma) %>% 
  append(ades.lm) %>% 
  append(ades.lme.c) %>% 
  append(ades.wi) %>% 
  append(ades.wi.c) %>% 
  append(ades.m)

addExperiments(algo.designs = ades, reg = reg)
submitJobs(x$job.id, resources=list(max.concurrent.jobs=2000))  
