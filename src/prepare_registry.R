# ##############################################################################
#
## Run all tests (vanilla use case) for a single simulation
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

sim.files <- list.files(paste0(temp.loc, 'simulations'), 
                        full.names = TRUE, recursive = TRUE)
sim.file <- sim.files[str_detect(sim.files, pattern = simulation)]
if (length(sim.file) != 1){
  stop("Cannot find the correct simulation!")
}
if (!file.exists(sim.file)){
  stop("No such simulation exists!")
}

sim.params <- h5read(sim.file, name='simulation_parameters/sim_params')
params <- yaml.load_file(here('create_simulations', "parameters.yaml"))

# tests <- list or vector
tests <- c('ANCOM', 'ANCOMBC', 'corncob', 'ZIBSeq', 'ZIBSeq-sqrt', 'DESeq2',
           'edgeR', 'metagenomeSeq', 'metagenomeSeq2', 'lm', 'wilcoxon',
           'limma', 'distinct', 'ZINQ', 'KS', 'ALDEx2')

# make experiment registry - might need to load pkgs here
if (!dir.exists(paste0(temp.loc, 'test_results_registries'))){
  dir.create(paste0(temp.loc, 'test_results_registires'))
}
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
temp <- h5ls(sim.file, recursive=FALSE)
groups <- setdiff(temp$name, c('original_data', 'simulation_parameters'))
# groups <- unique(str_remove(groups, '_rep[0-9]*$'))

### ADD PROBLEMS, ALGORITHMS, THEN EXPERIMENTS
sim.location <- sim.file

# add problems
addProblem(name = simulation,
           data = sim.location,
           fun = function(data, job) data,
           reg = reg)
  
# add algorithms (tests)
.f_apply <- function(data, job, instance, group, subsets, norm, test, ...){
  apply.test(sim.location = data, group=group, 
             subset=subsets, 
             norm=norm, 
             test=test)
}
walk(tests, ~ addAlgorithm(name = .,  fun = .f_apply))
  
# algorithm design
# a bit manual tweaking, but could work...
  
# slow tests
slow.tests <- c('ANCOM', 'corncob', 'ZIBSeq', 'ZIBSeq-sqrt', 
                'distinct', 'ALDEx2')
ades.s <- slow.tests %>%
  map(~ data.table::CJ(test = .,
                       group=groups,
                       norm = c('pass', 'rarefy'),
                       subsets=params$subset_size)) %>%
    set_names(slow.tests)
  
# very fast tests
v.fast.test <- c('limma', 'lm', 'KS')
ades.f <- v.fast.test %>% 
  map(~ data.table::CJ(test=.,
                       group=unique(str_remove(groups, '_rep[0-9]*$')),
                       norm = c('pass', 'TSS', 'TSS.log', 'clr',
                                'rclr', 'TSS.arcsin', 'rank', 'rarefy',
                                'rarefy.TSS', 'rarefy.TSS.log'),
                       subsets=NA_real_)) %>% 
  set_names(v.fast.test)
ades.w <- c('wilcoxon') %>% 
  map(~ data.table::CJ(test=.,
                       group=unique(str_remove(groups, '_rep[0-9]*$')),
                       norm = c('pass', 'TSS', 'clr', 'rclr', 'rarefy', 
                                'rarefy.TSS'),
                       subsets=NA_real_)) %>% 
  set_names(c('wilcoxon'))

  
# fast tests
rest <- setdiff(setdiff(tests, slow.tests), c(v.fast.test, 'wilcoxon'))
ades.r <- rest %>%
  map(~ data.table::CJ(test = .,
                       group=groups,
                       norm = c('pass', 'rarefy'),
                       subsets=NA_real_)) %>%
  set_names(rest)
  
ades <- append(ades.f, ades.r) %>% append(ades.s) %>% append(ades.w)
  
# add and submit jobs (no rep specified, default = 1)
addExperiments(algo.designs = ades, reg = reg)

