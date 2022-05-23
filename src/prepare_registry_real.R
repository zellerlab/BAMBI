# ##############################################################################
#
## Run tests for metacardis data with various drugs -- figure 4
#
# ##############################################################################

library(tidyverse)
library(batchtools)
library(here)

# load SIMBA
devtools::load_all(simba.loc)

args <- commandArgs(trailingOnly = TRUE)
folder <- args[1]
sim.pattern <- args[2]
test <- args[3]

if (!is.na(sim.pattern)){
  message("Taking simulation files with pattern: ", sim.pattern)
} else { 
  message("Taking all files in ", folder)
  list.files(here('simulations', folder)) }

if (!is.na(test)){
  message("Submitting jobs for the test: ", test)
} else { message("Submitting all different tests!") }

sim.files <- dir(here('simulations', folder), pattern = sim.pattern, 
                 full.names = TRUE)
# only works for 'pooled' sims currently
sim.names <- sim.files %>%
  str_split('_') %>%
  map_chr(~ magrittr::extract(., 5))

# HPC - jakob

# local - morgan
# job.registry <- here('test_results_registries',
#                      paste0(folder, '_test'))
# HPC - morgan
job.registry <- paste0('/scratch/essex/test_results_registries/',
                       folder, '_test')
message(job.registry)

if (!dir.exists(job.registry)) {
  reg <- makeExperimentRegistry(file.dir = job.registry,
                                work.dir = here(),
                                conf.file = here('cluster_config',
                                                 'batchtools_test_conf.R'),
                                make.default = TRUE)
} else {
  reg <- loadRegistry(file.dir = job.registry,
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

# add the simulations
for (sim.file in sim.files) {
  sim.name <- str_split(sim.file, '_') %>%
    map_chr(~ magrittr::extract(., 5)) %>%
    str_remove_all('.h5')
  addProblem(name = sim.name,
             data = here('simulations', folder, sim.file),
             fun = function(data, job) data,
             reg = reg) }

# updated algorithm to first check if confounder-aware test
.f_apply <- function(data, job, instance, group, type, 
                     subsets, norm, test, ...){
  if (grepl(test, pattern='_conf')){
    apply.test(sim.location = data,
               group=group,
               type=type,
               subset=subsets,
               norm=norm,
               test=gsub(x=test, pattern='_conf', replacement=''),
               conf='conf')
  } else {
    apply.test(sim.location = data,
               group=group,
               type=type,
               subset=subsets,
               norm=norm,
               test=test)
  }
}

# updated tests 
tests <- c('limma','limma_conf',
           'lm','lme_conf',
           'wilcoxon','wilcoxon_conf',
           'mdc-FE_conf')

walk(tests, ~ addAlgorithm(name = .,
                           fun = .f_apply))

# get metadata variables
temp <- rhdf5::h5ls(here('simulations',
                         folder,
                         sim.file), recursive = FALSE)
# params <- h5read(here('simulations', folder, sim.file),
#                  name = 'simulation_parameters')
# bias.terms <- c(0.5, 0.75)
# subsets <- c(50, 100, 200)
# groups <- setdiff(temp$name, c('original_data', 'simulation_parameters'))
# groups <- c('asa','clopidogrel','nitrate','dppiv','kspdiur','statine','antibiotics',
#            'beta_bl','ppi','metformin','antilipidtt','antithombo','heparine',
#            'insulin')
groups <- c('statine','antibiotics','beta_bl','ppi','metformin','asa','dppiv')

# add experiments
limma.tests <- c('limma','limma_conf')
ades.l <- limma.tests %>%
  map(~ data.table::CJ(test = .,
                       group = groups,
                       type = c('confounder','random'),
                       norm = c('pass','TSS.arcsin'),
                       subsets = NA_real_)) %>%
  set_names(limma.tests)
wcox.tests <- c('wilcoxon','wilcoxon_conf')
ades.w <- wcox.tests %>%
  map(~ data.table::CJ(test = .,
                       group = groups,
                       type = c('confounder','random'),
                       norm = c('pass','TSS'),
                       subsets = NA_real_)) %>%
  set_names(wcox.tests)

lm.tests <- c('lm','lme_conf','mdc-FE_conf') 
ades.lm <- lm.tests %>%
  map(~ data.table::CJ(test = .,
                       group = groups,
                       type = c('confounder','random'),
                       norm = c('pass','TSS.log'),
                       subsets = NA_real_)) %>%
  set_names(lm.tests)

ades <- append(ades.l, ades.w) %>% append(ades.lm)
ades <- ades[tests]

addExperiments(algo.designs = ades, reg = reg)
getJobPars() %>%
  unnest_wider(algo.pars) 
