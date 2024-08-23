# ##############################################################################
#
## Plot time and memory requirements for the different tools
#
# ##############################################################################

library("here")
library("tidyverse")
library("batchtools")
library("ggembl")
library("yaml")

reg <- here('test_results_registries', 'testing_time')
fn.time <- here('files', 'time_measurements_full.tsv')

if (file.exists(fn.time)){
  df.plot <- read_tsv(fn.time, col_types = cols())
} else if (dir.exists(reg)){
  loadRegistry(reg, writeable = TRUE)
  
  temp <- getJobTable() %>% 
    as_tibble()
  
  temp <- temp %>% 
    unnest_wider(algo.pars)
  
  df.missing <- temp %>% 
    filter(job.id %in% findNotDone()$job.id) %>% 
    group_by(algorithm) %>% 
    tally()
  if (nrow(df.missing) > 0){
    message("Some jobs did not finish yet!")
    for (a in seq_len(nrow(df.missing))){
      message("\t", df.missing$algorithm[a], ':', df.missing$n[a])
    }
  }
  df.plot <- temp %>% select(-c(prob.pars, resources))
  write_tsv(df.plot, fn.time)
} else if (!file.exists(fn.time)){
  params <- yaml.load_file(here('create_simulations', "parameters.yaml"))
  reg <- makeExperimentRegistry(file.dir = reg, 
                                work.dir = here(),
                                # conf.file = here('cluster_config',
                                                 # 'batchtools_test_conf.R'),
                                make.default = TRUE) 
  sim.file <- here('simulations', 'resampling', 'sim_Zeevi_WGS_all.h5')
  stopifnot(file.exists(sim.file))
  addProblem(name = 'time',
             data = sim.file,
             fun = function(data, job) data,
             reg = reg)
  tests <- c('ANCOM', 'ANCOM_old', 'ANCOMBC', 'corncob', 'ZIBSeq', 
             'metagenomeSeq', 'metagenomeSeq2', 'ZINQ',
             'DESeq2', 'edgeR', 'ALDEx2', 'distinct', 'limma', 'KS', 'lm',
             'wilcoxon', 'LinDA', 'fastANCOM', 'ZicoSeq', 't-test')
  # add algorithms (tests)
  .f_apply <- function(data, job, instance, group, subsets, norm, test, ...){
    apply.test(sim.location = data, group=group, 
               subset=subsets, 
               norm=norm, 
               test=test)
  }
  walk(tests, ~ addAlgorithm(name = ., 
                             fun = .f_apply))
  group <- 'ab4_prev2_rep10'
  ades <- tests %>%
    map(~ data.table::CJ(test = .,
                         group=group,
                         norm = c('pass'),
                         subsets=params$subset_size)) %>%
    set_names(tests)
  addExperiments(algo.designs = ades, reg = reg)
  submitJobs()
  # getJobTable() %>% 
  #   filter(job.id %in% findNotDone()$job.id) %>% 
  #   pull(algorithm) %>% table
  
  # once it is finished
  temp <- getJobTable() %>% 
    as_tibble()
  
  temp <- temp %>% 
    unnest_wider(algo.pars)
  
  df.missing <- temp %>% 
    filter(job.id %in% findNotDone()$job.id) %>% 
    group_by(algorithm) %>% 
    tally()
  if (nrow(df.missing) > 0){
    message("Some jobs did not finish yet!")
    for (a in seq_len(nrow(df.missing))){
      message("\t", df.missing$algorithm[a], ':', df.missing$n[a])
    }
  }
  df.plot <- temp %>% dplyr::select(-c(prob.pars, resources))
  write_tsv(df.plot, fn.time)
  
}

# ##############################################################################
# Plot time and memory

# group tests by something
# sort by time-requirements at 800 subset size
test.sort <- df.plot %>% 
  mutate(algorithm=str_remove(algorithm, '-sqrt')) %>% 
  mutate(algorithm=str_remove(algorithm, '_old')) %>% 
  filter(subsets==800) %>% 
  group_by(algorithm) %>% 
  slice(1) %>% 
  arrange(desc(time.running)) 
  
# test.colours <- viridisLite::viridis(n=nrow(test.sort))
# names(test.colours) <- test.sort$algorithm

test.colours <- unlist(read_yaml(here('files', 'test_colours.yml')))

df.plot <- df.plot %>% 
  mutate(type=case_when(str_detect(algorithm, '_old')~'old',
                        str_detect(algorithm, '-sqrt')~'sqrt',
                        TRUE~'standard')) %>% 
  mutate(algorithm_p=algorithm) %>% 
  mutate(algorithm=str_remove(algorithm, '-sqrt')) %>% 
  mutate(algorithm=str_remove(algorithm, '_old')) %>% 
  mutate(algorithm=factor(algorithm, levels = test.sort$algorithm))

g <- df.plot %>% 
  mutate(subsets=factor(subsets)) %>% 
  mutate(time_min=time.running/60) %>% 
  ggplot(aes(x=subsets, y=time_min, col=algorithm, linetype=type)) + 
    geom_line(aes(group=algorithm_p)) + 
    ylab("Run time for 50 repetitions [minutes]") + 
    xlab("Subset sample size") + 
    theme_publication(panel.grid = 'major') + 
    scale_colour_manual(values=test.colours, name='Algorithm') +
    scale_linetype_manual(values=c('standard'='solid', 'old'='dashed', 
                                   'sqrt'='dotted')) +
    NULL
ggsave(g, filename = here('figures', 'figure_vanilla', 'time_reqs.pdf'),
       width = 80, height = 80, units = 'mm', useDingbats=FALSE)

df.plot %>% 
  filter(subsets==800) %>% 
  select(algorithm_p, time.running) %>% 
  mutate(min=time.running/60) %>% 
  arrange(time.running)
