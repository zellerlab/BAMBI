# ##############################################################################
#
## Combine and plot the reality checks
#
# ##############################################################################

library("tidyverse")
library("batchtools")
library("devtools")
library("progress")
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
if (!file.exists(file.loc)){
  file.copy(from=sim.file, file.loc)
}

sim.params <- h5read(file.loc, name='simulation_parameters/sim_params')

# get groups
temp <- h5ls(file.loc, recursive = FALSE)
groups <- setdiff(temp$name, c('original_data', 'simulation_parameters'))

# load the registry
job.registry <- paste0(temp.loc, 'reality_checks_registries/', simulation)
if (!dir.exists(job.registry)){
  stop("Job registry does not exist!")
}
loadRegistry(file.dir=job.registry)

# check that all jobs are done
stats <- getStatus()
if (stats$defined == 0){
  stop("Jobs have not been defined yet!")
}
if (stats$defined != stats$done){
  stop("Some jobs haven't finished yet!")
}

# intialize results folder
res.dir <- paste0(temp.loc, 'reality_checks_results/', simulation)
if (!dir.exists(here(res.dir))){
  dir.create(here(res.dir))
}

# ##############################################################################
# load all results
all.jobs <- getJobPars() %>%
  as_tibble() %>%
  unnest(cols=job.pars) %>%
  unnest(cols=job.pars) %>%
  mutate(ab=str_extract(job.pars, 'ab[0-9]+')) %>%
  mutate(ab=as.numeric(str_remove(ab, 'ab')))
if ('prev.scale' %in% names(sim.params)){
  all.jobs <- all.jobs %>%
    mutate(prev=str_extract(job.pars, 'prev[0-9]+')) %>%
    mutate(prev=as.numeric(str_remove(prev, 'prev')))
} else {
  all.jobs <- all.jobs %>%
    mutate(prev=0)
}

df.auc <- tibble(ab=integer(0), prev=integer(0), rep=integer(0),
                 AUC=double(0), type=character(0), F_value=double(0), 
                 R_value=double(0), P_value=double(0))
df.sparsity <- list()
df.variance <- list()
df.correlation <- tibble(ab=integer(0), prev=integer(0), rep=integer(0),
                         corr=double(0))
pb <- progress_bar$new(total=nrow(all.jobs) * sim.params$repeats)
for (a in seq_len(nrow(all.jobs))){
  info <- all.jobs %>%
    select(-job.pars) %>%
    slice(a)

  temp <- loadResult(info$job.id)

  for (i in seq_along(temp)){
	res <- temp[[i]]
    # AUC
    df.auc <- df.auc %>%
      bind_rows(tibble(ab=info$ab, prev=info$prev, rep=i,
                       AUC=as.numeric(res$ML), type='ML')) %>% 
      bind_rows(tibble(ab=info$ab, prev=info$prev, rep=i,
                       AUC=as.numeric(res$Bray_Curtis$AUC), 
                       type='Bray-Curtis', F_value=res$Bray_Curtis$F,
                       R_value=res$Bray_Curtis$R2, 
                       P_value=res$Bray_Curtis$P)) %>% 
      bind_rows(tibble(ab=info$ab, prev=info$prev, rep=i,
                       AUC=as.numeric(res$Euclidean$AUC), 
                       type='Euclidean', F_value=res$Euclidean$F,
                       R_value=res$Euclidean$R2, 
                       P_value=res$Euclidean$P))
      
    # correlation
    df.correlation <- df.correlation %>% 
      add_row(ab=info$ab, prev=info$prev, rep=i, corr=res$Correlation)
    
    # sparsity
    df.sparsity[[(length(df.sparsity) + 1)]] <- res$sparsity %>%
      select(-group) %>%
      mutate(ab=info$ab, prev=info$prev, rep=i)
    # variance
    df.variance[[(length(df.variance) + 1)]] <- res$variance %>%
      mutate(ab=info$ab, prev=info$prev, rep=i)  
    pb$tick()
  }
}

df.auc <- df.auc %>%
  mutate(simulation=simulation)
df.sparsity <- bind_rows(df.sparsity) %>%
  mutate(simulation=simulation)
df.variance <- bind_rows(df.variance) %>%
  mutate(simulation=simulation)

# get measures for the original dataset
original.data <- reality.check(file.loc, 'original')
df.sparsity <- df.sparsity %>%
  bind_rows(original.data$sparsity %>%
              select(-group) %>%
              mutate(ab=0, prev=0, rep=0) %>%
              mutate(simulation=simulation))
df.variance <- df.variance %>%
  bind_rows(original.data$variance %>%
              mutate(ab=0, prev=0, rep=0) %>%
              mutate(simulation=simulation) %>%
              mutate(fc.sim=NA_real_, prev.shift.sim=NA_real_))
df.correlation <- df.correlation %>% 
  bind_rows(tibble(ab=0, prev=0, rep=0, corr=original.data$Correlation))

# save results
write_tsv(df.auc, file = here(res.dir,'auc_all.tsv'))
write_tsv(df.correlation, file = here(res.dir,'correlation.tsv'))
write_tsv(df.sparsity, file = here(res.dir, 'sparsity.tsv'))
write_tsv(df.variance, file = here(res.dir, 'variance.tsv'))

# ##############################################################################
# actually plot a few PCoAs?

pdf(here(res.dir, 'pco_plots.pdf'),
    useDingbats = FALSE, width = 6, height = 5)
if ('prev.scale' %in% names(sim.params)){
  for (a in seq_along(sim.params$ab.scale)){
    for (b in seq_along(sim.params$prev.scale)){
      for (dist in c('bray', 'log-euclidean')){
        r <- sample(sim.params$repeats, 1)
        x <- pcoa.plot(sim.location = file.loc,
                       group=paste0('ab', a, '_prev', b, '_rep', r),
                       distance = dist)
        x <- x + ggtitle(paste0(dist, "-distance, Ab: ",
                                sim.params$ab.scale[a], ' Prev: ',
                                sim.params$prev.scale[b]))
        print(x)
      }
    }
  }

} else {
  for (a in seq_along(sim.params$ab.scale)){
    for (dist in c('bray', 'log-euclidean')){
      r <- sample(sim.params$repeats, 1)
      x <- pcoa.plot(sim.location = file.loc,
                     group=paste0('ab', a, '_rep', r),
                     distance = dist)
      x <- x + ggtitle(paste0(dist, "-distance, Ab: ", sim.params$ab.scale[a]))
      print(x)
    }
  }
}
dev.off()

# ##############################################################################
# Plot the rest

# AUCs
auc.colours <- c('#D0DEBB', '#6CC24A', '#009F4D', '#007A53', '#115740')
g.auc <- df.auc %>%
  group_by(ab, prev, type) %>%
  summarise(x=mean(AUC), xs=sd(AUC), .groups = 'drop') %>%
  ggplot(aes(x=as.factor(ab), y=as.factor(prev), fill=x)) +
    geom_tile() +
    theme_minimal() +
    scale_fill_gradientn(colours = auc.colours, name='AUC',
                         limits=c(0.48, 1)) +
    geom_text(aes(label=paste0(sprintf(fmt='%.2f', x), "\u00B1",
                               sprintf(fmt='%.2f', xs))),
              col='white') +
    xlab('Abundance scaling') +
    ylab('Prevalence scaling') +
    facet_grid(type~., scales = 'free', space = 'free') +
    theme(panel.grid = element_blank())
ggsave(g.auc, filename = here(res.dir, 'auc_plot.pdf'),
       useDingbats=FALSE, width = 9,
       height = ifelse('prev.scale' %in% names(sim.params), 10, 4))

# sparsity
sparsity.colours <- c('#8BB8E8', "#307FE2", "#003DA5")
n <- h5read(file=file.loc, name='/original_data/filt_features') %>% nrow
g.sparsity <- df.sparsity %>%
  mutate(L0=L0/n) %>%
  pivot_longer(cols=c(Gini, L0),
               names_to = 'type', values_to = 'values') %>%
  group_by(ab, prev, type) %>%
  summarise(x=mean(values), xs=sd(values), .groups = 'drop') %>%
  ggplot(aes(x=as.factor(ab), y=as.factor(prev), fill=x)) +
    geom_tile() +
    theme_minimal() +
    scale_fill_gradientn(colours = sparsity.colours, name='Sparsity fraction',
                         limits=c(0, 1)) +
    geom_text(aes(label=paste0(sprintf(fmt='%.2f', x), "\u00B1",
                               sprintf(fmt='%.2f', xs))),
              col='white') +
    xlab('Abundance scaling') +
    ylab('Prevalence scaling') +
    facet_grid(type~., scales = 'free', space = 'free') +
    theme(panel.grid = element_blank())
ggsave(g.sparsity, filename = here(res.dir, 'sparsity_plot.pdf'),
       useDingbats=FALSE, width = 9,
       height = ifelse('prev.scale' %in% names(sim.params), 6, 3))

# variance
var.sim <- df.variance %>%
  filter(ab!=0) %>%
  arrange(ab, prev, rep, var.rel)
var.original <- df.variance %>%
  filter(ab==0) %>%
  filter(features %in% var.sim$features) %>%
  arrange(var.rel)
df.plot <- var.sim %>%
  mutate(var.rel.original=
           rep(var.original$var.rel,
               max(var.sim$ab)*max(var.sim$rep)*max(max(var.sim$prev), 1))) %>%
  mutate(diff=log10(var.rel) - log10(var.rel.original))

# median difference in variance by group
colour.bar <- c("#A6093D", "#E40046", "#E58F9E", "white",
                "#8BB8E8", "#307FE2", "#003DA5")
g <- df.plot %>%
  group_by(selection, ab, prev, feat.type) %>%
  summarise(x=median(diff), .groups = 'drop') %>%
  mutate(label=case_when(abs(x) > 0.5~sprintf(fmt='%.2f', x),
                         TRUE~'')) %>%
  mutate(feat.type=factor(feat.type, levels = c('low', 'middle', 'high'))) %>%
  ggplot(aes(x=as.factor(ab), y=as.factor(prev), fill=x))  +
    geom_tile() +
    facet_grid(selection~feat.type, scales = 'free', space = 'free') +
    scale_fill_gradientn(colours=colour.bar,
                         limits=c(-3, 3),
                         name='Median change\nin variance') +
    theme_minimal() +
    geom_text(aes(label=label), col='white') +
    xlab('Abundance scaling') +
    ylab('Prevalence scaling') +
    theme(panel.grid=element_blank())
ggsave(g, filename = here(res.dir, 'variance_by_group_and_type.pdf'),
       useDingbats=FALSE, width = 11,
       height = ifelse('prev.scale' %in% names(sim.params), 12, 5))
g <- df.plot %>%
  group_by(selection, ab, prev) %>%
  summarise(x=median(diff), .groups='drop') %>%
  mutate(label=case_when(abs(x) > 0.5~sprintf(fmt='%.2f', x),
                         TRUE~'')) %>%
  ggplot(aes(x=as.factor(ab), y=as.factor(prev), fill=x))  +
    geom_tile() +
    facet_grid(.~selection, scales = 'free', space = 'free') +
    scale_fill_gradientn(limits=c(-2.5, 2.5),
                         colours = colour.bar,
                         name='Median change\nin variance') +
    theme_minimal() + theme(panel.grid=element_blank()) +
    geom_text(aes(label=label), col='white') +
    xlab('Abundance scaling') +
    ylab('Prevalence scaling')
ggsave(g, filename = here(res.dir, 'variance_by_group.pdf'),
       useDingbats=FALSE, width = 11,
       height = ifelse('prev.scale' %in% names(sim.params), 9, 3))

# scatter plots
g <- df.plot %>%
  group_by(ab, prev, selection, features) %>%
  summarise(var.rel.sim=median(var.rel),
            var.rel.original=median(var.rel.original), .groups='drop') %>%
  ggplot(aes(x=log10(var.rel.original), y=log10(var.rel.sim),
             colour=selection)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(alpha=0.3) +
    theme_bw() +
    xlab('Original feature variance') +
    ylab('Simulated feature variance')  +
    scale_colour_manual(values=c('#707372', '#E40046')) +
    theme(panel.grid.minor = element_blank()) +
    facet_grid(prev~ab)
ggsave(g, filename = here(res.dir, 'variance_scatter_plots.pdf'),
       useDingbats=FALSE, width = 9,
       height = ifelse('prev.scale' %in% names(sim.params), 10, 4))

# effect sizes
g.effect <- df.variance %>%
  filter(ab != 0) %>%
  mutate(prev.shift.sim=
           case_when(is.infinite(prev.shift.sim)~
                       max(prev.shift.sim[is.finite(prev.shift.sim)]),
                     TRUE~prev.shift.sim)) %>%
  group_by(ab, prev) %>%
  filter(rep==sample(max(rep))) %>%
  ungroup() %>%
  ggplot(aes(x=abs(fc.sim), y=abs(log2(prev.shift.sim)),
             col=selection, alpha=selection)) +
    geom_point() +
    facet_grid(prev~ab) +
    xlab('Absolute fold change') +
    ylab('log2(Prevalence shift)') +
    theme_bw() +
    scale_colour_manual(values=c('#707372', "#307FE2")) +
    scale_alpha_manual(values=c(0.2, 0.6)) +
    theme(panel.grid.minor=element_blank())
ggsave(g.effect, filename = here(res.dir, 'effect_size.pdf'),
       useDingbats=FALSE, width = 9,
       height = ifelse('prev.scale' %in% names(sim.params), 10, 4))
