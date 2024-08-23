# ##############################################################################
#
## Figure 2: Reality assessment for the vanilla benchmark
#
# ##############################################################################

library("here")
library("tidyverse")
library("vroom")
library("ggembl")
library("ggrepel")

embl.greens <- function(x){
  greens <- colorRamp(colors=c('#D0DEBB', '#6CC24A', '#18974C', 
                               '#007B53', '#0A5032'))
  rgb(greens(seq(from=0, to=1, length.out=x)), maxColorValue = 255)
}
embl.greys <- function(x){
  greens <- colorRamp(colors=c('#D0D0CE', '#A8A99E', '#707372', 
                               '#54585A', '#373A36'))
  rgb(greens(seq(from=0, to=1, length.out=x)), maxColorValue = 255)
}

sim.colours <- c('reference'='#707372', 'MMH'='#A6093D50', 'Weiss'='#A6093D',
                 'betabin'='#E58F9E', 'negbin'='#D41645', 'nbc'='#D4164550',
                 'dir'='#FFD020', 'SimMSeq'='', 'sparseDOSSA'='#FDB71A', 
                 'all'='#3B6FB6')
sim.colours <- c('reference'='#707372', 'MMH'='#651C3250', 'Weiss'='#651C32',
                 'betabin'='#E58F9E', 'negbin'='#8C1515', 'nbc'='#8C151575',
                 'dir'='#E04F39', 'SimMSeq'='#FFBF00', 'sparseDOSSA'='#E98300', 
                 'all'='#3B6FB6')

enframe(sim.colours) %>%
  mutate(name=factor(name, levels = name)) %>% 
  ggplot(aes(x=name, y=1, fill=name)) + 
    geom_tile() + 
    scale_fill_manual(values = sim.colours) + 
    theme_void()

# order all other simulations for the Supplement
sim.levels <- c('reference', 
                'MMH', 'Weiss', 'betabin', 'negbin', 'nbc', 'dir', 
                'SimMSeq', 'sparseDOSSA',
                'all', 'abundance', 'high', 'middle', 
                'low_middle', 'low', 'inverse-abundance',
                'compositional')

# ##############################################################################
# combine all measures
comb <- FALSE
if (comb){
  all.sims <- list.files(here('reality_checks'), pattern = 'sim_')
  extra <- c('abundance', 'high', 'middle',  'low_middle', 'low', 
             'inverse-abundance', 'compositional')
  extra.sims <- all.sims[str_detect(all.sims, paste(extra, collapse = '|'))]
  all.sims <- setdiff(all.sims, extra.sims)
  message("+ Found ", length(all.sims),  ' and ', length(extra.sims), 
          " additional different simulations with reality assessment!")
  # AUC
  message('++ Reading AUC tables')
  auc.all <- map(all.sims, .f = function(x){
    temp <- read_tsv(here('reality_checks', x, 'auc_all.tsv'), 
                     col_types = cols())
    info <- unique(temp$simulation)
    if (length(info) < 1){stop("Check AUC file for: ", x)}
    info <- str_split(info, pattern = '_')[[1]]
    temp <- temp %>% 
      select(-simulation) %>% 
      mutate(dataset=info[2], assay=info[3], sim_type=info[4])}) %>% 
    bind_rows()
  write_tsv(auc.all, here('reality_checks', 'auc_all.tsv'))
  # extra stuff
  auc.extra <- map(extra.sims, .f = function(x){
    temp <- read_tsv(here('reality_checks', x, 'auc_all.tsv'), 
                     col_types = cols())
    info <- unique(temp$simulation)
    if (length(info) < 1){stop("Check AUC file for: ", x)}
    info <- str_split(info, pattern = '_')[[1]]
    temp <- temp %>% 
      select(-simulation) %>% 
      mutate(dataset=info[2], assay=info[3], sim_type=info[4])}) %>% 
    bind_rows()
  write_tsv(auc.extra, here('reality_checks', 'auc_extra.tsv'))
  
  # Variance
  message('++ Reading variance tables')
  var.all <- map(all.sims, .f = function(x){
    temp <- read_tsv(here('reality_checks', x, 'variance.tsv'), 
                     col_types = cols())
    info <- unique(temp$simulation)
    if (length(info) < 1){stop("Check variance file for: ", x)}
    info <- str_split(info, pattern = '_')[[1]]
    temp <- temp %>% 
      select(-simulation) %>% 
      mutate(dataset=info[2], assay=info[3], sim_type=info[4])}) %>% 
    bind_rows()
  reference <- var.all %>% 
    group_by(dataset, assay) %>% 
    filter(ab==0) %>% 
    filter(sim_type=='MMH') %>% 
    mutate(sim_type='reference')
  var.all <- var.all %>% 
    filter(ab!=0)
  var.all <- bind_rows(var.all, reference)
  write_tsv(var.all, here('reality_checks', 'variance_all.tsv'))
  var.extra <- map(extra.sims, .f = function(x){
    temp <- read_tsv(here('reality_checks', x, 'variance.tsv'), 
                     col_types = cols())
    info <- unique(temp$simulation)
    if (length(info) < 1){stop("Check variance file for: ", x)}
    info <- str_split(info, pattern = '_')[[1]]
    temp <- temp %>% 
      select(-simulation) %>% 
      mutate(dataset=info[2], assay=info[3], sim_type=info[4])}) %>% 
    bind_rows()
  reference <- var.extra %>% 
    group_by(dataset, assay) %>% 
    filter(ab==0) %>% 
    filter(sim_type=='abundance') %>% 
    mutate(sim_type='reference')
  var.extra <- var.extra %>% 
    filter(ab!=0)
  var.extra <- bind_rows(var.extra, reference)
  write_tsv(var.extra, here('reality_checks', 'variance_extra.tsv'))
  
  # Sparsity
  message('++ Reading sparsity tables')
  sparsity.all <- map(all.sims, .f = function(x){
    temp <- read_tsv(here('reality_checks', x, 'sparsity.tsv'), 
             col_types = cols())
    info <- unique(temp$simulation)
    if (length(info) < 1){stop("Check sparsity file for: ", x)}
    info <- str_split(info, pattern = '_')[[1]]
    temp <- temp %>% 
      select(-simulation) %>% 
      mutate(dataset=info[2], assay=info[3], sim_type=info[4])
    }) %>% 
    bind_rows() 
  reference <- sparsity.all %>% 
    group_by(dataset, assay) %>% 
    filter(ab==0) %>% 
    filter(sim_type=='MMH') %>% 
    mutate(sim_type='reference')
  sparsity.all <- sparsity.all %>% 
    filter(ab!=0)
  sparsity.all <- bind_rows(sparsity.all, reference)
  write_tsv(sparsity.all, here('reality_checks', 'sparsity_all.tsv'))
  sparsity.extra <- map(extra.sims, .f = function(x){
    temp <- read_tsv(here('reality_checks', x, 'sparsity.tsv'), 
                     col_types = cols())
    info <- unique(temp$simulation)
    if (length(info) < 1){stop("Check sparsity file for: ", x)}
    info <- str_split(info, pattern = '_')[[1]]
    temp <- temp %>% 
      select(-simulation) %>% 
      mutate(dataset=info[2], assay=info[3], sim_type=info[4])
  }) %>% 
    bind_rows() 
  reference <- sparsity.extra %>% 
    group_by(dataset, assay) %>% 
    filter(ab==0) %>% 
    filter(sim_type=='MMH') %>% 
    mutate(sim_type='reference')
  sparsity.extra <- sparsity.extra %>% 
    filter(ab!=0)
  sparsity.extra <- bind_rows(sparsity.extra, reference)
  write_tsv(sparsity.extra, here('reality_checks', 'sparsity_extra.tsv'))
  
  # correlation
  message('++ Reading correlation tables')
  corr.all <- map(all.sims, .f=function(x){
    fn <- here('reality_checks', x, 'correlation.tsv')
    info <- str_split(x, pattern = '_')[[1]]
    if (file.exists(fn)){
      read_tsv(fn, col_types = cols()) %>% 
        mutate(dataset=info[2], assay=info[3], sim_type=info[4])} 
    else {NULL}}) %>% 
    bind_rows()
  reference <- corr.all %>% 
    group_by(dataset, assay) %>% 
    filter(ab==0) %>% 
    filter(sim_type=='MMH') %>% 
    mutate(sim_type='reference')
  corr.all <- corr.all %>% 
    filter(ab!=0)
  corr.all <- bind_rows(corr.all, reference)
  write_tsv(corr.all, here('reality_checks', 'correlation_all.tsv'))
}

# ##############################################################################
# Compare implantation vs other sims
auc.all <- read_tsv(here('reality_checks', 'auc_all.tsv'), col_types = cols())
sparsity.all <- read_tsv(here('reality_checks', 'sparsity_all.tsv'), 
                         col_types = cols())
variance.all <- read_tsv(here('reality_checks', 'variance_all.tsv'), 
                         col_types = cols())
# plot for Figure 1b
g.var <- variance.all %>% 
  filter(dataset=='Zeevi', assay=='WGS') %>% 
  filter(prev %in% c(0, 1)) %>% 
  filter(ab %in% c(0, 4)) %>% 
  filter(rep %in% c(0, 2)) %>% 
  mutate(sim_type=factor(sim_type, levels = sim.levels)) %>% 
  ggplot(aes(x=log10(var.rel), col=sim_type)) + 
    geom_density()+
    scale_colour_manual(values=sim.colours) + 
    theme_publication() + 
    xlab('log10(Feature variance)') + 
    ylab("Density")
ggsave(g.var, filename = here('figures', 'figure_reality', 
                              'var_density_main.pdf'),
       width = 5, height = 3, useDingbats=FALSE)  
# corresponding supplement
g.var.all <- variance.all %>% 
  filter(assay!='KEGG') %>% 
  mutate(d=paste0(dataset, '_', assay)) %>% 
  filter(prev %in% c(0, 1)) %>% 
  filter(ab %in% c(0, 4)) %>% 
  filter(rep %in% c(0, 2)) %>% 
  mutate(sim_type=factor(sim_type, levels = sim.levels)) %>% 
  ggplot(aes(x=log10(var.rel), col=sim_type)) + 
    geom_density()+
    scale_colour_manual(values=sim.colours) + 
    theme_publication() + 
    xlab('log10(Feature variance)') + 
    ylab("Density") + 
    facet_wrap(d~., scales = 'free')
ggsave(g.var.all, filename = here('figures', 'figure_reality', 
                                  'var_density_all.pdf'),
       width = 7, height = 5, useDingbats=FALSE)  

g.var.all.ridges <- variance.all %>% 
  filter(assay!='KEGG') %>% 
  mutate(d=paste0(dataset, '_', assay)) %>% 
  filter(prev %in% c(0, 1)) %>% 
  filter(ab %in% c(0, 4)) %>% 
  filter(rep %in% c(0, 2)) %>% 
  filter(var.rel!=0) %>% 
  mutate(sim_type=case_when(sim_type=='dirmult'~'dir', TRUE~sim_type)) %>% 
  mutate(sim_type=factor(sim_type, levels = sim.levels)) %>%
  ggplot(aes(x=log10(var.rel), y=sim_type, fill=sim_type)) + 
    ggridges::geom_density_ridges() +
    scale_fill_manual(values=sim.colours, guide='none') + 
    theme_bw() + 
    xlab('log10(Feature variance)') + ylab('') + 
    theme(panel.grid = element_blank(), axis.text.y = element_blank(), 
          axis.ticks.y=element_blank()) + 
    facet_wrap(~d, scales = 'free_x')
ggsave(g.var.all.ridges, filename = here('figures', 'figure_reality', 
                                  'var_density_all_ridges.pdf'),
       width = 7, height = 5, useDingbats=FALSE)  
  
# how about the overdispersion stuff?
impl.feat <- variance.all %>% 
  filter(dataset=='Zeevi', assay=='WGS', sim_type=='all') %>% 
  filter(selection) %>% 
  arrange(mean.rel) %>% 
  mutate(mean.bin=cut(log10(mean.rel), 50))
tmp <- variance.all %>% 
  filter(dataset=='Zeevi', assay=='WGS', sim_type=='all') %>% 
  filter(!selection) %>% 
  mutate(mean.bin=cut(log10(mean.rel), breaks=levels(impl.feat$mean.bin)))


# plot for Figure 1c
# combine AUC/PERMANOVA/Sparsity?
# adjust sparsity
sparsity.prop <- sparsity.all %>% 
  mutate(dataset=paste0(dataset, '-', assay)) %>% 
  left_join(variance.all %>% 
              filter(sim_type=='reference') %>% 
              mutate(dataset=paste0(dataset, '-', assay)) %>% 
              group_by(dataset) %>% tally, by="dataset") %>% 
  mutate(value=L0/n)
g.measures <- auc.all %>% 
  filter(ab < 6) %>% 
  filter(dataset=='Zeevi', assay=='WGS') %>% 
  filter(type %in% c("ML", 'Euclidean')) %>% 
  mutate(AUC=case_when(AUC<0.5~0.5, TRUE~AUC)) %>% 
  mutate(value=case_when(type=='ML'~AUC, 
                         type=='Euclidean'~log10(F_value+1))) %>% 
  select(sim_type, value, type) %>% 
  bind_rows(sparsity.prop %>% 
              filter(ab<6) %>% 
              filter(dataset=='Zeevi-WGS') %>% 
              transmute(value, sim_type, type='Sparsity')) %>% 
  mutate(sim_type=factor(sim_type, levels = sim.levels)) %>% 
  ggplot(aes(x=sim_type, y=value, col=type)) + 
    geom_boxplot() + 
    facet_grid(type~., scales = 'free') + 
    theme_publication() + 
    xlab('') + ylab('') + 
    scale_colour_manual(values=c('#BE5400', '#00B7BD', '#ABD037'))
ggsave(g.measures, filename = here('figures', 'figure_reality', 
                                   'measures.pdf'),
       width = 5, height = 4, useDingbats=FALSE)

# corresponding supplement
g.measures.all <- auc.all %>% 
  filter(ab < 6) %>% 
  filter(assay!='KEGG') %>% 
  filter(type %in% c("ML", 'Euclidean')) %>% 
  mutate(AUC=case_when(AUC<0.5~0.5, TRUE~AUC)) %>% 
  mutate(value=case_when(type=='ML'~AUC, 
                         type=='Euclidean'~log10(F_value+1))) %>% 
  mutate(dataset=paste0(dataset, '-', assay)) %>% 
  select(sim_type, value, type, dataset) %>% 
  bind_rows(sparsity.prop %>% 
              filter(ab<6) %>% 
              filter(assay!='KEGG') %>% 
              transmute(value, sim_type, type='Sparsity', dataset)) %>% 
  mutate(sim_type=factor(sim_type, levels = sim.levels)) %>% 
  ggplot(aes(x=sim_type, y=value, col=type)) + 
  geom_boxplot() + 
  facet_grid(type~dataset, scales = 'free_y') + 
  theme_publication() + 
  xlab('') + ylab('') + 
  scale_colour_manual(values=c('#BE5400', '#00B7BD', '#ABD037')) + 
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(g.measures.all, filename = here('figures', 'figure_reality', 
                                       'measures_all.pdf'),
       width = 9, height = 4, useDingbats=FALSE)


# supplement: mean-var relationship
# mean-variance relationship
ref <- variance.all %>% filter(sim_type=='reference') %>% 
  mutate(dataset=paste0(dataset, '-', assay)) %>% 
  filter(assay!='KEGG')
ref %>% 
  filter(dataset=='Zeevi-WGS') %>% 
  ggplot(aes(x=log10(mean.rel), y=log10(var.rel))) + 
  stat_density2d_filled(breaks=c(0, 0.001, 0.05, 0.1,
                                 0.15,0.2, 0.25, 0.3, 3)) +
  # facet_grid(~dataset) +
  scale_fill_manual(values=c('#FFFFFF00', embl.greys(12)[1:7]),
                    guide=FALSE) +
  theme_presentation(panel.grid='major') + 
  geom_point(alpha=0.6) + 
  xlab('log10(mean relative abundance)') + 
  ylab('log10(variance)')

g <- variance.all %>% 
  mutate(dataset=paste0(dataset, '-', assay)) %>% 
  filter(assay!='KEGG') %>% 
  filter(prev%in%c(0,1), ab%in%c(0,4), rep%in%c(0,2)) %>% 
  mutate(sim_type=factor(sim_type, levels = sim.levels)) %>% 
  arrange(selection) %>% 
  ggplot(aes(x=log10(mean.rel), y=log10(var.rel), col=selection)) + 
  stat_density2d_filled(data=ref %>% select(-sim_type), colour='#FFFFFF00',
                        breaks=c(0, 0.001, 0.05, 0.1, 0.15,0.2, 0.25, 0.3, 3)) +  
  geom_point(alpha=0.3) +
  facet_grid(dataset+assay~sim_type) + 
  scale_colour_manual(values=c('#54585A', '#D41645')) +
  scale_fill_manual(values=c('#FFFFFF00', embl.greys(12)[1:7])) + 
  theme_publication(panel.grid = 'major')
ggsave(g, filename = here('figures','figure_reality', 
                          'mean_var_relationship.pdf'),
       useDingbats=FALSE, width = 380, height = 300, units='mm')

# mean-variance for more effect sizes
g.mean.var.ef <- variance.all %>% 
  filter(dataset=='Zeevi', assay=='WGS', sim_type=='all') %>% 
  filter(rep==100, prev==2) %>% 
  # filter(ab%in%c(2, 3, 4, 5), prev==2) %>%  
  arrange(selection) %>%
  ggplot(aes(x=log10(mean.rel), y=log10(var.rel), col=selection)) + 
  facet_grid(prev~ab) +
  stat_density2d_filled(data=ref %>% 
                          filter(dataset=='Zeevi-WGS', assay=='WGS') %>% 
                          select(-sim_type, -ab, -prev), colour='#FFFFFF00',
                        breaks=c(0, 0.001, 0.05, 0.1, 0.15,0.2, 0.25, 0.3)) +  
  geom_point(alpha=0.3) +
  scale_colour_manual(values=c('#54585A', '#D41645')) +
  scale_fill_manual(values=c('#FFFFFF00', embl.greys(12)[1:6])) + 
  theme_publication(panel.grid = 'major')
ggsave(g.mean.var.ef, width = 8, height = 3, useDingbats=FALSE,
       filename = here('figures', 'misc', 
                       'mean_variance_effect_sizes.pdf'))

## quantify this somehow!
# calculate the dispersion for each group (real data vs simulated data)
if (!file.exists(here('reality_checks','permdisp.tsv'))){
  temp <- variance.all %>% 
    filter(assay!='KEGG') %>% 
    mutate(dataset=paste0(dataset, '-', assay)) %>% 
    mutate(group=paste0('ab', ab, '_prev', prev, '_rep', rep)) %>% 
    group_by(dataset, sim_type, group) %>% 
    mutate(var.rel=log10(var.rel+1e-12), mean.rel=log10(mean.rel+1e-7)) %>% 
    mutate(centroid.dist=sqrt((var.rel-mean(var.rel))^2 + 
                                (mean.rel-mean(mean.rel))^2)) %>% 
    mutate(group=case_when(group=='ab0_prev0_rep0'~'reference', TRUE~group))

  df.test.permdisp <- tibble(group=character(0), pval=double(0), 
                             f_value=double(0), simulation=character(0),
                             dataset=character(0))
  for (d in unique(temp$dataset)){
    message(d)
    red <- temp %>% filter(dataset==d)
    for (sim.type in unique(red$sim_type)){
      message(sim.type)
      red.sim <- red %>% 
        filter(dataset==d, sim_type %in% c(sim.type, 'reference'))
      pb <- progress::progress_bar$new(total=length(unique(red.sim$group)) -1)
      for (x in setdiff(unique(red.sim$group), 'reference')){
        red.group <- red.sim %>% filter(group %in% c(x, 'reference'))
        fit <- lm(red.group$centroid.dist~red.group$group)
        res <- anova(fit)
        df.test.permdisp <- df.test.permdisp %>% 
          add_row(group=x, pval=res[1,5], f_value=res[1,4],
                  simulation=sim.type, dataset=d)
        pb$tick()
      }
    }
  }

  write_tsv(df.test.permdisp, file = here('reality_checks','permdisp.tsv'))
} else {
  df.test.permdisp <- read_tsv(here('reality_checks','permdisp.tsv'))
}

df.test.permdisp %>% 
  mutate(group=str_remove(group, '_rep[0-9]*')) %>% 
  group_by(dataset, simulation, group) %>% 
  summarise(m=mean(pval < 0.05), .groups = 'drop') %>% 
  ggplot(aes(x=group, y=simulation, fill=m)) + 
    geom_tile() +
    facet_wrap(~dataset, scales='free_x') + 
    theme_publication() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    xlab('') + ylab('') + 
    ggthemes::scale_fill_gradient_tableau()


# supplement: effect sizes across extra simulations
# all effect sizes for the supplement
var.extra <- read_tsv(here('reality_checks', 'variance_extra.tsv'), 
                      col_types = cols())
var.sum <- variance.all %>% 
  filter(dataset=='Zeevi', assay=='WGS') %>% 
  filter(sim_type!='reference') %>% 
  bind_rows(var.extra %>% filter(sim_type!='reference')) %>% 
  group_by(sim_type, ab, prev, rep, selection) %>% 
  summarise(fc=mean(abs(fc.sim), na.rm=TRUE), 
            prev.shift=mean(abs(prev.shift.sim), na.rm=TRUE))

# viridis colour-scheme?
colour.scheme <- c(
  "#FDE725FF", "#fef5a7", "#fdf07c", "#fdeb50", "#FDE725FF",
  "#8FD744FF", "#d2efb4", "#bbe78e", "#a5df69", "#8FD744FF",
  "#35B779FF", "#aee2c9", "#85d3ae", "#5dc593", "#35B779FF",
  "#21908CFF", "#a6d2d1", "#79bcba", "#4da6a3", "#21908CFF",
  "#31688EFF", "#acc2d1", "#83a4bb", "#5a86a4", "#31688EFF",
  "#443A83FF", "#b4b0cd", "#8e88b4", "#69619b", "#443a83FF",
  "#440154FF", "#b499ba", "#8e6698", "#693376", "#440154FF"
)

g <- var.sum %>% 
  mutate(group=paste0(ab, '-', prev)) %>% 
  mutate(sim_type=factor(sim_type, sim.levels)) %>% 
  ggplot(aes(x=fc, y=prev.shift, col=group, shape=selection)) + 
  geom_point() + 
  facet_wrap(~sim_type) + 
  theme_publication() + 
  xlab("Absolute fold change") + 
  ylab("Absolute prevalence shift") + 
  scale_colour_manual(values = colour.scheme, 
                      labels = c('1-0', '2-0', '3-0', '4-0', 
                                 '5-0', '6-0', '7-0'),
                      breaks = c('1-0', '2-0', '3-0', '4-0', 
                                 '5-0', '6-0', '7-0')) + 
  scale_shape_manual(values=c(3, 16))
ggsave(g, filename = here('figures', 'figure_reality', 'all_effect_size.pdf'),
       useDingbats=FALSE, units = 'mm', height = 180, width = 180)

# how about background features with higher fold change?
tmp <- variance.all %>% 
  filter(dataset=='Zeevi', sim_type=='all', assay=='WGS')
lvls <- tmp %>% 
  group_by(rep, ab, prev) %>% 
  filter(selection) %>% 
  mutate(fc.sim=abs(fc.sim)) %>% 
  summarise(lvl=quantile(fc.sim, probs=0.15), .groups='drop')
tmp %>% 
  full_join(lvls) %>% 
  mutate(higher=abs(fc.sim) > lvl) %>% 
  group_by(ab, prev, rep, selection) %>% 
  summarise(prob=sum(higher)/n()) %>% 
  filter(!selection) %>% 
  ggplot(aes(x=as.factor(ab), y=prob, fill=as.factor(prev))) + 
    geom_boxplot() + 
    theme_embl(panel.grid='major') + 
    xlab('Abundance scaling effect size') + 
    ylab('Fraction of features with gFC > 0.2 quantile of implanted features')
tmp %>% 
  full_join(lvls) %>% 
  mutate(higher=abs(fc.sim) > lvl) %>% 
  group_by(ab, prev, rep, selection, higher) %>% 
  summarise(n=n(), .groups='drop') %>% 
  mutate(type=paste0(selection, higher)) %>% 
  filter(type %in% c('FALSETRUE', 'TRUETRUE')) %>% 
  select(-c(higher, type)) %>% 
  pivot_wider(names_from = selection, values_from = n) %>% 
  mutate(ratio=`FALSE`/`TRUE`) %>% 
  filter(ab==5) %>% 
  ggplot(aes(x=as.factor(ab), y=log10(ratio), fill=as.factor(prev))) + 
    geom_boxplot() + 
    theme_embl(panel.grid='major') + 
    xlab('Abundance scaling effect size') + 
    ylab('Ratio')

tmp.auc <- tmp %>% 
  group_by(ab, prev, rep) %>% 
  summarise(auc=as.numeric(pROC::auc(response=selection, 
                                     predictor=abs(fc.sim))))
g <- tmp.auc %>% 
  # mutate(auc=case_when(auc < 0.5~0.5, TRUE~auc)) %>% 
  ggplot(aes(x=as.factor(ab), y=auc, fill=as.factor(prev))) + 
    geom_boxplot() + 
    theme_publication(panel.grid='major') + 
    xlab('Abundance scaling effect size') + 
    ylab('AUC') +
    scale_x_discrete(labels=c('1.0', '1.25', '1.5', '2.0', '5.0', 
                              '10.0', '20.0')) +
    scale_fill_manual(values=c("#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF"), 
                      name='Prevalence', labels=c('0.0', '0.1', '0.2', '0.3'))
ggsave(g, filename=here('figures', 'figure_reality', 'gfc_auc.pdf'), 
                        width = 6, height=3)

# ##############################################################################
# additional revisions: SPSimSeq






# ##############################################################################
# example features for main figure

# create measures for real datasets
if (!file.exists(here('files', 'real_effect_sizes.tsv'))){
  # CRC
  library("coin")
  feat.crc <- read.table(here('data', 'motus_crc_rel_meta.tsv'), sep='\t',
                         stringsAsFactors = FALSE, check.names = FALSE,
                         quote = '', comment.char = '') %>% 
    as.matrix()
  meta.crc <- read_tsv(here('data', 'meta_crc.tsv'), col_types = cols())
  res.crc <- vapply(rownames(feat.crc), FUN=function(x){
    case <- feat.crc[x,meta.crc %>% filter(Group=='CRC') %>% pull(Sample_ID)]
    ctr <- feat.crc[x,meta.crc %>% filter(Group=='CTR') %>% pull(Sample_ID)]
    q.case <- quantile(log10(case+1e-05), probs = seq(0.05, 0.95, by=0.05))
    q.ctr <- quantile(log10(ctr+1e-05), probs = seq(0.05, 0.95, by=0.05))
    wilcox.test <- wilcox_test(feat.crc[x,meta.crc$Sample_ID]~
                                 as.factor(meta.crc$Group)|
                                 as.factor(meta.crc$Study))
    p.val <- pvalue(wilcox.test)
    return(c('fc'=mean(q.case-q.ctr), 
             'prev.diff'=mean(case != 0) - mean(ctr != 0),
             'p.val'=p.val))
  }, FUN.VALUE=double(3)) %>% t() %>% as_tibble(rownames='species') %>% 
    mutate(disease='CRC') %>% 
    mutate(p.adj=p.adjust(p.val, method='fdr'))
  # CD
  feat.cd <- read.table(here('data', 'motus_cd_rel_meta.tsv'), sep='\t',
                        stringsAsFactors = FALSE, check.names = FALSE,
                        quote = '', comment.char = '') %>% 
    as.matrix()
  meta.cd <- read_tsv(here('data', 'meta_cd.tsv'), col_types = cols())
  res.cd <- vapply(rownames(feat.cd), FUN=function(x){
    case <- feat.cd[x,meta.cd %>% filter(Group=='CD') %>% pull(Sample_ID)]
    ctr <- feat.cd[x,meta.cd %>% filter(Group=='CTR') %>% pull(Sample_ID)]
    q.case <- quantile(log10(case+1e-05), probs = seq(0.05, 0.95, by=0.05))
    q.ctr <- quantile(log10(ctr+1e-05), probs = seq(0.05, 0.95, by=0.05))
    wilcox.test <- wilcox_test(feat.cd[x,meta.cd$Sample_ID]~
                                 as.factor(meta.cd$Group)|
                                 as.factor(meta.cd$Study))
    p.val <- pvalue(wilcox.test)
    return(c('fc'=mean(q.case-q.ctr), 
             'prev.diff'=mean(case != 0) - mean(ctr != 0),
             'p.val'=p.val))
  }, FUN.VALUE=double(3)) %>% t() %>% as_tibble(rownames='species') %>% 
    mutate(disease='CD') %>% 
    mutate(p.adj=p.adjust(p.val, method='fdr'))
  real.effect.size <- bind_rows(res.crc, res.cd)
  write_tsv(real.effect.size, file = here('files', 'real_effect_sizes.tsv'))
} else {
  real.effect.size <- read_tsv(here('files', 'real_effect_sizes.tsv'))
  feat.crc <- read.table(here('data', 'motus_crc_rel_meta.tsv'), sep='\t',
                         stringsAsFactors = FALSE, check.names = FALSE,
                         quote = '', comment.char = '') %>% 
    as.matrix()
  meta.crc <- read_tsv(here('data', 'meta_crc.tsv'), col_types = cols())
  feat.cd <- read.table(here('data', 'motus_cd_rel_meta.tsv'), sep='\t',
                        stringsAsFactors = FALSE, check.names = FALSE,
                        quote = '', comment.char = '') %>% 
    as.matrix()
  meta.cd <- read_tsv(here('data', 'meta_cd.tsv'), col_types = cols())
}

# example measures
sim.measures <- variance.all %>% 
  filter(sim_type=='all') %>% 
  filter(dataset=='Zeevi', assay=='WGS') %>% 
  bind_rows(var.extra %>% 
              filter(sim_type=='low.h')) %>% 
  filter(ab==4, prev==2, rep==86)

real.effect.size %>% 
  ggplot(aes(x=abs(fc), y=abs(prev.diff), col=disease)) + 
    geom_point() + 
    xlab('Absolute generalized fold change') + 
    ylab('Absolute prevalence difference') + 
    theme_bw() +
    facet_grid(~disease) 

sim.measures %>% 
  ggplot(aes(x=abs(fc.sim), y=abs(prev.shift.sim), col=selection)) + 
    geom_point() + 
    xlab('Absolute generalized fold change') + 
    ylab('Absolute prevalence shift)') + 
    theme_bw() +
    facet_grid(~sim_type) + 
    scale_colour_manual(values=c('#E4004675', '#FFA30075'))


real.effect.size <- real.effect.size %>% 
  mutate(colour_group=case_when(p.adj < 1e-20~'real-3',
                                p.adj < 1e-10~'real-2',
                                p.adj < 1e-05~'real-1',
                                TRUE~'real'))
# plot CRC
real.effect.size %>% 
  filter(disease=='CRC') %>% 
  mutate(species=str_remove(species, 
                            '\\[(ref|meta)_mOTU_v25_[0-9]{5}\\]$')) %>% 
  mutate(highlight=case_when(colour_group=='real-3'~species,
                             str_detect(species, 'nucleatum')~species,
                             TRUE~'')) %>% 
  ggplot(aes(x=abs(fc), y=abs(prev.diff), col=colour_group)) + 
    geom_point() + 
    xlab('Absolute generalized fold change') + 
    ylab('Absolute prevalence difference)') + 
    theme_bw() + 
    geom_text_repel(aes(label=highlight), min.segment.length = 0) + 
    scale_colour_manual(values=c('#70737250', '#8BB8E8', 
                                 '#307FE2', '#003DA5'))

# plot CD
real.effect.size %>% 
  filter(disease=='CD') %>% 
  mutate(species=str_remove(species, 
                            '\\[(ref|meta)_mOTU_v25_[0-9]{5}\\]$')) %>% 
  mutate(fc=abs(fc), prev.diff=abs(prev.diff)) %>% 
  mutate(highlight=case_when(
    fc > 0.75 & prev.diff < 0.2 ~ species,
    prev.diff > 0.5~species,
    fc > 1.4 ~ species,
    TRUE~'')) %>% 
  mutate(unc=str_detect(species, 
                        'Clostridiales species incertae sedis')) %>% 
  ggplot(aes(x=fc, y=prev.diff, col=colour_group)) + 
    geom_point(aes(shape=unc)) + 
    xlab('Absolute generalized fold change') + 
    ylab('Absolute prevalence difference)') + 
    theme_bw() + 
    geom_text_repel(aes(label=highlight), min.segment.length = 0) + 
    scale_colour_manual(values=c('#70737250', '#8BB8E8', 
                                 '#307FE2', '#003DA5')) + 
    scale_shape_manual(values=c(16, 1))

# select species
# CRC
sp.1 <- 'Parvimonas micra [ref_mOTU_v25_04287]'
sp.2 <- 'Fusobacterium nucleatum subsp. nucleatum [ref_mOTU_v25_01003]'
# IBD
sp.3 <- 'Clostridiales species incertae sedis [meta_mOTU_v25_12635]'
sp.4 <- '[Ruminococcus] gnavus [ref_mOTU_v25_01594]'#'Proteobacteria sp. [ref_mOTU_v25_00095]'

# plot together
sim.measures %>% 
  mutate(fc=abs(fc.sim)) %>% 
  mutate(prev.diff=abs(prev.shift.sim)) %>% 
  transmute(prev.diff, fc, colour_group=case_when(selection~'implanted', 
                                                   TRUE~'background'), 
            disease=sim_type, type='simulated', 
            selection=as.character(selection)) %>% 
  bind_rows(real.effect.size %>% mutate(fc=abs(fc), 
                                    prev.diff=abs(prev.diff),
                                    selection='real')) %>% 
  arrange(p.adj, selection) %>% 
  mutate(highlight=species %in% c(sp.1, sp.2, sp.3, sp.4)) %>% 
  ggplot(aes(x=fc, y=prev.diff, fill=colour_group)) + 
    geom_point(aes(colour=highlight), shape=21) + 
    facet_grid(~disease) + 
    xlab('Abs. gFC') + 
    ylab("Abs log2(prevalence shift)") + 
    theme_publication() + 
    scale_fill_manual(values=c('#A8A99E50', '#A1BE1F75', '#70737250', 
                               '#6CC24A75', '#18974C75', '#0A503275')) + 
  scale_colour_manual(values=c('#FFFFFF00', 'black'))


real.effect.size %>% filter(species %in% c(sp.1, sp.2, sp.3, sp.4))
sim.measures %>% filter(sim_type=='all') %>% 
  filter(abs(fc.sim) > 0.5) %>% arrange(abs(prev.shift.sim))
sim.measures %>% filter(sim_type=='low.h') %>% 
  arrange(desc(abs(fc.sim)))
sim.measures %>% filter(sim_type=='low.h') %>% 
  filter(abs(fc.sim) < 0.25) %>% arrange(desc(abs(prev.shift.sim)))

# all
sp.5 <- 'bact_otu_13944'
sp.6 <- 'bact_otu_4000'
# low
sp.7 <- 'bact_otu_4145'
sp.8 <- 'bact_otu_350'#'bact_otu_350'

sim.measures %>% 
  mutate(highlight=features %in% c(sp.5, sp.6, sp.7, sp.8)) %>% 
  ggplot(aes(x=abs(fc.sim), y=abs(prev.shift.sim), col=selection)) + 
  geom_point()  +
  xlab('Absolute generalized fold change') + 
  ylab('Absolute prevalence shift)') + 
  theme_bw() +
  facet_grid(~sim_type) + 
  scale_colour_manual(values=c('#E4004675', '#FFA30075')) +
  geom_point(data=. %>% filter(highlight), col='black')

g <- sim.measures %>% 
  mutate(fc=abs(fc.sim)) %>% 
  mutate(prev.diff=abs(prev.shift.sim)) %>% 
  transmute(prev.diff, fc, colour_group=case_when(selection~'implanted', 
                                                  TRUE~'background'), 
            disease=sim_type, type='simulated', species=features,
            selection=as.character(selection)) %>% 
  bind_rows(real.effect.size %>% mutate(fc=abs(fc), 
                                        prev.diff=abs(prev.diff),
                                        selection='real')) %>% 
  arrange(p.adj, selection) %>% 
  mutate(highlight=species %in% c(sp.1, sp.2, sp.3, sp.4,
                                  sp.5, sp.6, sp.7, sp.8)) %>% 
  ggplot(aes(x=fc, y=prev.diff, fill=colour_group)) + 
  geom_point(aes(colour=highlight), shape=21) + 
  facet_grid(~disease) + 
  xlab('Abs. gFC') + 
  ylab("Abs log2(prevalence shift)") + 
  theme_publication() + 
  scale_fill_manual(values=c('#A8A99E50', '#A1BE1F75', '#70737250', 
                             '#6CC24A75', '#18974C75', '#0A503275')) + 
  scale_colour_manual(values=c('#FFFFFF00', 'black'))
ggsave(g, filename = here('figures', 'figure_reality', 
                          'markers_comparison.pdf'),
       width = 7, height = 4, useDingbats=FALSE)
# show some examples
# get original data
# get sim_all features
temp <- rhdf5::h5read('./simulations/resampling/sim_Zeevi_WGS_all.h5', 
                      name='ab4_prev2_rep86')
sim.feat.all <- temp$features
rownames(sim.feat.all) <- temp$feature_names
colnames(sim.feat.all) <- temp$sample_names
sim.feat.all <- prop.table(sim.feat.all, 2)
meta.all <- tibble(Sample_ID=temp$sample_names, 
                   Group=temp$labels, type='all')

temp <- rhdf5::h5read('./simulations/resampling/sim_Zeevi_WGS_low.h5',
                      name='ab4_prev2_rep86')
sim.feat.low <- temp$features
rownames(sim.feat.low) <- temp$feature_names
colnames(sim.feat.low) <- temp$sample_names
sim.feat.low <- prop.table(sim.feat.low, 2)
meta.low <- tibble(Sample_ID=temp$sample_names, 
                   Group=-temp$labels, type='low.h5')

df.plot <-
  t(feat.crc[c(sp.1, sp.2),meta.crc$Sample_ID]) %>% 
  as_tibble(rownames='Sample_ID') %>% 
  pivot_longer(-Sample_ID, names_to = 'species', values_to = 'rel.ab') %>% 
  left_join(meta.crc %>% mutate(type='CRC'), by='Sample_ID') %>% 
  bind_rows(
    t(feat.cd[c(sp.3, sp.4),meta.cd$Sample_ID]) %>% 
      as_tibble(rownames='Sample_ID') %>% 
      pivot_longer(-Sample_ID, names_to = 'species', values_to = 'rel.ab') %>% 
      left_join(meta.cd %>% mutate(type='CD'), by='Sample_ID')) %>% 
  bind_rows(
    t(sim.feat.all[c(sp.5, sp.6),meta.all$Sample_ID]) %>%
      as_tibble(rownames='Sample_ID') %>%
      pivot_longer(-Sample_ID, names_to = 'species', values_to = 'rel.ab') %>%
      left_join(meta.all %>% 
                  mutate(Group=as.character(Group)) %>% 
                  mutate(Group=case_when(Group==1~'Group1', TRUE~'Group2')), 
                by='Sample_ID') %>% 
      mutate(Group=case_when(species==sp.5 & Group=='Group1'~'Group2', 
                             species==sp.5 & Group=='Group2'~'Group1', 
                             TRUE~Group))) %>%
  bind_rows(
    t(sim.feat.low[c(sp.7, sp.8),meta.low$Sample_ID]) %>%
      as_tibble(rownames='Sample_ID') %>%
      pivot_longer(-Sample_ID, names_to = 'species', values_to = 'rel.ab') %>%
      left_join(meta.low %>%
                  mutate(Group=as.character(Group)) %>%
                  mutate(Group=case_when(Group==1~'Group1', TRUE~'Group2')),
                by='Sample_ID') %>%
      mutate(Group=case_when(species==sp.7 & Group=='Group1'~'Group2',
                             species==sp.7 & Group=='Group2'~'Group1',
                             TRUE~Group)))

g <- df.plot %>% 
  mutate(rel.ab=log10(rel.ab+1e-05)) %>% 
  mutate(Group=factor(Group, levels = c('CTR', 'CD', 'CRC', 
                                        'Group1', 'Group2'))) %>%
  mutate(type=factor(type, levels = c('CRC', 'low.h5', 'CD', 'all'))) %>% 
  ggplot(aes(x=species, y=rel.ab, fill=Group)) + 
    geom_quantile_box(alpha=0.75) + 
    geom_jitter(aes(col=Group), 
                position = position_jitterdodge(jitter.width = 0.2)) + 
    scale_fill_manual(values=c('#70737250', '#18974C','#18974C', 
                               '#F49E17', '#734595')) + 
    scale_colour_manual(values=c('#707372', '#18974C','#18974C', 
                                 '#F49E17', '#734595')) + 
    theme_publication() + 
    xlab('') + 
    theme(axis.ticks.x = element_blank()) + 
    facet_grid(~type, scales = 'free')
ggsave(g, filename = here('figures', 'figure_reality', 'selected_markers.pdf'),
       useDingbats=FALSE, units = 'mm', height = 70, width = 180)
