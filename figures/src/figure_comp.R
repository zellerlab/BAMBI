# ##############################################################################
#
## Supplementary Figure: Reality assessment and vanilla results for 
##    compositionally-confounded simulations
#
# ##############################################################################

library("tidyverse")
library("here")
library("ggembl")
library("ggrepel")

# ##############################################################################
# AUC
auc.comp <- read_tsv(here('reality_checks', 'compositional', 'auc_all.tsv'))
auc.rest <- read_tsv(here('reality_checks', 'auc_all.tsv'))

g <- auc.comp %>% 
  bind_rows(auc.rest) %>% 
  filter(type=='ML') %>% 
  filter(simulation %in% c('sim_compositional', 'sim_all', 
                           'sim_MMH', 'sim_negbin')) %>% 
  ggplot(aes(x=simulation, y=AUC, fill=as.factor(ab))) + 
  geom_boxplot() +
  theme_publication() + 
  ylab("log10(F value)") +
  scale_fill_manual(values=viridis::viridis(n=7), name='Abundance')


auc.values.rest <- c(1, 1.25, 1.5, 2, 5, 10, 20)
auc.values.comp <- c(1, 1.25, 1.5, 2, 2.5, 3, 3.5, 4, 5, 10, 20)
g <- auc.comp %>% 
  mutate(ab.value=auc.values.comp[ab]) %>% 
  bind_rows(auc.rest %>% 
              mutate(ab.value=auc.values.rest[ab])) %>% 
  filter(type=='Euclidean') %>% 
  filter(simulation %in% c('compositional', 'sim_all', 
                           'sim_MMH', 'sim_negbin')) %>% 
  mutate(prev=paste0('prev_', prev)) %>% 
  ggplot(aes(x=simulation, y=log10(F_value + 1), fill=as.factor(ab.value))) + 
    geom_boxplot() +
    theme_publication() + 
    ylab("log10(F value)") +
    scale_fill_manual(values=viridis::viridis(n=11), name='Abundance') +
    NULL
ggsave(g, filename = here('figures', 'figure_comp', 'auc.pdf'),
       width = 7, height = 4, useDingbats=FALSE)

# ##############################################################################
# effect size plot

effect.size.all <- read_tsv(here('reality_checks', 'sim_all', 'variance.tsv'))
effect.size.comp <- read_tsv(here('reality_checks', 
                                  'compositional', 'variance.tsv'))

g <- effect.size.all %>% 
  filter(prev==2) %>% 
  filter(rep==48) %>% 
  bind_rows(effect.size.comp %>% filter(prev==2) %>% filter(rep==48)) %>% 
  ggplot(aes(x=as.factor(ab), y=fc.sim, fill=selection)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.08),
                size=0.7) + 
    facet_grid(~simulation, scales = 'free_x', space = 'free') +
    theme_publication(panel.grid = 'major_y') + 
    xlab('Abundance scaling') + 
    ylab('Generalized fold change') + 
    scale_fill_manual(values=alpha(c('#707372', '#3B6FB6'), alpha=0.7), 
                      name='Implanted')
ggsave(g, filename = here('figures', 'figure_comp', 'effect_size.pdf'),
       width = 7, height = 4, useDingbats=FALSE)


# ##############################################################################
# results from vanilla testing

# take a specific group/sample size
fn.final <- here('figures', 'figure_comp', 'test_results.tsv')
if (!file.exists(fn.test.results)){
  fn.test.results <- here('test_results', 'compositional_all.tsv')
  df.res.all <- read_tsv(fn.test.results) %>% 
    mutate(group=str_remove(group, '_rep[0-9]*$')) %>%
    mutate(subset=str_remove(subset, 'subset_')) %>% 
    group_by(test, group, subset) %>% 
    summarise(mean.FDR=mean(FDR), sd.FDR=sd(FDR),
              mean.R=mean(R), sd.R=sd(R),
              mean.auroc=mean(auroc), sd.auroc=sd(auroc)) %>% 
    separate(group, into=c('ab', 'prev'), sep='_')
  write_tsv(df.res.all, file=fn.final)
} else {
  df.res.all <- read_tsv(fn.final)
}

test.colours <- yaml::yaml.load_file(here('files', 'test_colours.yml'))

g <- df.res.all %>% 
  mutate(ab=factor(ab, levels = paste0('ab', seq_len(11)))) %>% 
  filter(prev=='prev2') %>% 
  select(test, ab, mean.FDR, mean.R, mean.auroc, subset) %>% 
  mutate(mean.PR=1-mean.FDR) %>% 
  pivot_longer(cols=c(mean.FDR, mean.R, mean.auroc, mean.PR)) %>% 
  filter(name!='mean.FDR') %>% 
  mutate(lt=test=='metagenomeSeq2') %>% 
  ggplot(aes(x=ab, y=value)) + 
    geom_line(aes(group=test, col=test, lty=lt))  +
    facet_grid(subset~name) +
    xlab('Abundance scaling') + 
    ylab('') +
    scale_x_discrete(labels=c('1', '1.25', '1.5', '2', '2.5', '3', '3.5', 
                              '4', '5', '10', '20')) +
    theme_publication(panel.grid = 'major_y') + 
    scale_colour_manual(values = unlist(test.colours))
    
ggsave(g, filename = here('figures', 'figure_comp', 'all_measures.pdf'),
       width = 5, height = 4, useDingbats=FALSE)



df.res.all %>% 
  mutate(ab=factor(ab, levels = paste0('ab', seq_len(11)))) %>% 
  arrange(ab) %>% 
  filter(prev=='prev2') %>% 
  select(test, ab, mean.FDR, mean.R, mean.auroc, subset) %>% 
  filter(subset==400) %>% 
  mutate(mean.PR=1-mean.FDR) %>% 
  ggplot(aes(x=mean.R, y=mean.PR, col=test)) + 
    geom_point(aes(size=ab)) + 
    geom_line(aes(group=test)) + 
    scale_colour_manual(values = unlist(test.colours))
