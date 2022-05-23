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

# decide on the included simulations for the main figure
sims <- c('CRC'='CRC', 'IBD'='IBD',
          'reference'='reference',
          'sim_MMH'='multinomial',
          'sim_negbin'='negative binonmial', 
          'sim_betabin'='beta binonmial', 
          'sim_dirichlet'='dirichlet', 
          'sim_sparseDOSSA'='sparseDOSSA', 
          'sim_all'='resampling')

# order all other simulations for the Supplement
sim.levels <- c('reference', 
                'sim_MMH', 'sim_Weiss', 'sim_betabin', 'sim_negbin', 
                'sim_negbin_cor', 'sim_dirichlet', 'sim_sparseDOSSA', 
                'sim_all', 'sim_abundance', 'sim_high', 'sim_middle', 
                'sim_low_middle', 'sim_low', 'sim_inverse-abundance',
                'sim_compositional')

# ##############################################################################
# combine all measures
comb <- FALSE
if (comb){
  all.sims <- list.files(here('reality_checks'), pattern = 'sim_')
  # AUC
  auc.all <- map(all.sims, .f = function(x){
    read_tsv(here('reality_checks', x, 'auc_all.tsv'))}) %>% 
    bind_rows()
  write_tsv(auc.all, here('reality_checks', 'auc_all.tsv'))
  
  # Variance
  var.all <- map(all.sims, .f = function(x){
    read_tsv(here('reality_checks', x, 'variance.tsv'))}) %>% 
    bind_rows()
  reference <- var.all %>% 
    filter(ab==0) %>% 
    filter(simulation=='sim_MMH') %>% 
    mutate(simulation='reference')
  var.all <- var.all %>% 
    filter(ab!=0)
  var.all <- bind_rows(var.all, reference)
  write_tsv(var.all, here('reality_checks', 'variance_all.tsv'))
  
  # Sparsity
  sparsity.all <- map(all.sims, .f = function(x){
    read_tsv(here('reality_checks', x, 'sparsity.tsv'))}) %>% 
    bind_rows()
  reference <- sparsity.all %>% 
    filter(ab==0) %>% 
    filter(simulation=='sim_MMH') %>% 
    mutate(simulation='reference')
  sparsity.all <- sparsity.all %>% 
    filter(ab!=0)
  sparsity.all <- bind_rows(sparsity.all, reference)
  write_tsv(sparsity.all, here('reality_checks', 'sparsity_all.tsv'))
  
  # correlation
  corr.all <- map(all.sims, .f=function(x){
    fn <- here('reality_checks', x, 'correlation.tsv')
    if (file.exists(fn)){
      read_tsv(fn) %>% 
        mutate(simulation=x)} 
    else {NULL}}) %>% 
    bind_rows()
  reference <- corr.all %>% 
    filter(ab==0) %>% 
    filter(simulation=='sim_MMH') %>% 
    mutate(simulation='reference')
  corr.all <- corr.all %>% 
    filter(ab!=0)
  corr.all <- bind_rows(corr.all, reference)
  write_tsv(corr.all, here('reality_checks', 'correlation_all.tsv'))
}

# ##############################################################################
# AUC
auc.all <- read_tsv(here('reality_checks', 'auc_all.tsv'))

# PERMANOVA plot
g <- auc.all %>% 
  filter(!is.na(F_value)) %>% 
  filter(simulation %in% names(sims)) %>% 
  mutate(sim=sims[simulation]) %>% 
  mutate(sim=factor(sim, levels = sims)) %>% 
  filter(type=='Euclidean') %>% 
  # filter(ab==4) %>% 
  ggplot(aes(x=sim, y=log10(F_value+1))) + 
    geom_boxplot() + 
    xlab('') + 
    theme_publication()
ggsave(g, filename = here('figures', 'figure_reality', 'permanova.pdf'),
       useDingbats=FALSE, width = 61, height = 30, units = 'mm')
# all simulations
g <- auc.all %>% 
  mutate(simulation=factor(simulation, levels = sim.levels)) %>% 
  filter(!is.na(F_value)) %>% 
  filter(prev < 2) %>% 
  ggplot(aes(x=simulation, y=log10(F_value+1), fill=as.factor(ab))) + 
    geom_boxplot() + 
    xlab('') + ylab('log10(F value)') +
    theme_publication() + 
    facet_grid(type~., scales = 'free') + 
    scale_fill_manual(values=c("#D0DEBB", "#90CC71", "#54B64B", "#009F4D", 
                               "#038651", "#086E4D", "#115740"),
                      name='Effect size') + 
    theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(g, filename = here('figures', 'figure_reality', 'permanova_all.pdf'),
       useDingbats=FALSE, width = 160, height = 90, units = 'mm')

# AUC plot
g <- auc.all %>% 
  filter(type=='ML') %>% 
  filter(simulation %in% names(sims)) %>% 
  mutate(sim=sims[simulation]) %>% 
  mutate(sim=factor(sim, levels = sims)) %>% 
  mutate(AUC=case_when(AUC <= 0.5~1-AUC, TRUE~AUC)) %>% 
  # filter(ab==4) %>% 
  # filter(prev %in%  c(0, 3)) %>%
  ggplot(aes(x=sim, y=AUC)) + 
    geom_boxplot() + 
    xlab('') + 
    theme_publication()
ggsave(g, filename = here('figures', 'figure_reality', 'auc.pdf'),
       useDingbats=FALSE, width = 61, height = 30, units = 'mm')
# all AUC
g <- auc.all %>% 
  mutate(simulation=factor(simulation, levels = sim.levels)) %>% 
  filter(type=='ML') %>% 
  filter(prev < 2) %>%
  mutate(AUC=case_when(AUC <= 0.5 ~ 1-AUC, TRUE~AUC)) %>% 
  ggplot(aes(x=simulation, y=AUC, fill=as.factor(ab))) + 
    geom_boxplot() + 
    xlab('') + 
    theme_publication() + 
    facet_grid(type~., scales = 'free') + 
    scale_fill_manual(values=c("#D0DEBB", "#90CC71", "#54B64B", "#009F4D", 
                               "#038651", "#086E4D", "#115740"),
                      name='Effect size') + 
    theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(g, filename = here('figures', 'figure_reality', 'auc_all.pdf'),
       useDingbats=FALSE, width = 160, height = 70, units = 'mm')

# ##############################################################################
# Sparsity plot
sparsity.all <- read_tsv(here('reality_checks', 'sparsity_all.tsv'))

# selected
g <- sparsity.all %>% 
  filter(simulation %in% names(sims)) %>% 
  mutate(sim=sims[simulation]) %>% 
  mutate(sim=factor(sim, levels = c(sims))) %>% 
  ggplot(aes(x=sim, y=L0)) +   
    geom_boxplot() + 
    theme_publication() + 
    xlab('') + 
    ylab('Sparsity (L0)')
ggsave(g, filename = here('figures', 'figure_reality', 'sparsity.pdf'),
       useDingbats=FALSE, width = 120, height = 50, units = 'mm')

# all 
g <- sparsity.all %>% 
  filter(prev < 2) %>% 
  mutate(simulation=factor(simulation, levels = sim.levels)) %>% 
  ggplot(aes(x=simulation, y=L0, fill=as.factor(ab))) + 
    geom_boxplot() + 
    xlab('') + 
    theme_publication() + 
    scale_fill_manual(values=c("grey", 
                               "#D0DEBB", "#90CC71", "#54B64B", "#009F4D", 
                               "#038651", "#086E4D", "#115740"),
                      name='Effect size') + 
    theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(g, filename = here('figures', 'figure_reality', 'sparsity_all.pdf'),
       useDingbats=FALSE, width = 160, height = 70, units = 'mm')

# ##############################################################################
# Correlation plot
correlation.all <- read_tsv(here('reality_checks', 'correlation_all.tsv'))

# selected
g <- correlation.all %>% 
  filter(simulation %in% names(sims)) %>% 
  mutate(sim=sims[simulation]) %>% 
  mutate(sim=factor(sim, levels = c(sims))) %>% 
  ggplot(aes(x=sim, y=corr)) +   
    geom_boxplot() + 
    theme_publication() + 
    xlab('') + 
    ylab('Feature correlation')
ggsave(g, filename = here('figures', 'figure_reality', 'correlation.pdf'),
       useDingbats=FALSE, width = 120, height = 50, units = 'mm')

# all 
g <- correlation.all %>% 
  filter(prev < 2) %>% 
  mutate(simulation=factor(simulation, levels = sim.levels)) %>% 
  ggplot(aes(x=simulation, y=corr, fill=as.factor(ab))) + 
  geom_boxplot() + 
  xlab('') + 
  theme_publication() + 
  scale_fill_manual(values=c("grey", 
                             "#D0DEBB", "#90CC71", "#54B64B", "#009F4D", 
                             "#038651", "#086E4D", "#115740"),
                    name='Effect size') + 
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(g, filename = here('figures', 'figure_reality', 'correlation_all.pdf'),
       useDingbats=FALSE, width = 160, height = 70, units = 'mm')

# ##############################################################################
# Variance plot
var.all <- read_tsv(here('reality_checks', 'variance_all.tsv'))

g <- var.all %>% 
  filter(simulation %in% names(sims)) %>% 
  filter(prev %in% c(0, 1)) %>% 
  filter(ab %in% c(0, 4)) %>% 
  filter(rep %in% c(0, 2)) %>% 
  mutate(sim=sims[simulation]) %>% 
  mutate(sim=factor(sim, levels = c(sims))) %>% 
  mutate(var.rel=log10(var.rel)) %>% 
  ggplot(aes(x=var.rel, col=sim)) + 
    geom_density() + 
    theme_publication() + 
    xlab('log10(Feature variance)')
ggsave(g, filename = here('figures', 'figure_reality', 'var_density.pdf'), 
       useDingbats=FALSE, width = 120, height = 50, units = 'mm')

# supplement: all feature variance things
g <- var.all %>% 
  filter(prev < 2) %>% 
  mutate(simulation=factor(simulation, levels = sim.levels)) %>% 
  mutate(var.rel=log10(var.rel)) %>% 
  ggplot(aes(x=as.factor(ab), y=var.rel, fill=selection)) + 
    geom_boxplot() + 
    xlab('Abundance scaling effect size') +
    ylab('log(Feature variance)') + 
    theme_publication() + 
    facet_wrap(~simulation) +
    scale_fill_manual(values=c('#70737275', '#E4004675'))
ggsave(g, filename = here('figures', 'figure_reality', 'var_all.pdf'), 
       useDingbats=FALSE, width = 180, height = 180, units = 'mm')

# ##############################################################################
# Effect size plot
ext.measures <- read_tsv('./files/external_measures.tsv')
# var.all <- vroom::vroom(here('reality_checks', 'variance_all.tsv'))

# example measures
sim.measures <- var.all %>% 
  filter(simulation%in%c('sim_all', 'sim_low')) %>% 
  filter(ab==4, prev==2, rep==22)
# rm(var.all)

ext.measures %>% 
  ggplot(aes(x=abs(fc), y=abs(prev.shift), col=group)) + 
    geom_point() + 
    xlab('Absolute generalized fold change') + 
    ylab('Absolute prevalence shift)') + 
    theme_bw() +
    facet_grid(~group) + 
    scale_colour_manual(values=c('#E4004675', '#FFA30075'))

sim.measures %>% 
  ggplot(aes(x=abs(fc.sim), y=abs(prev.shift.sim), col=selection)) + 
    geom_point() + 
    xlab('Absolute generalized fold change') + 
    ylab('Absolute prevalence shift)') + 
    theme_bw() +
    facet_grid(~simulation) + 
    scale_colour_manual(values=c('#E4004675', '#FFA30075'))


ext.measures <- ext.measures %>% 
  mutate(colour_group=case_when(p.adj < 1e-20~'real-3',
                                p.adj < 1e-10~'real-2',
                                p.adj < 1e-05~'real-1',
                                TRUE~'real'))
# plot CRC
ext.measures %>% 
  filter(group=='CRC') %>% 
  mutate(species=str_remove(species, 
                            '\\[(ref|meta)_mOTU_v25_[0-9]{5}\\]$')) %>% 
  mutate(highlight=case_when(colour_group=='real-3'~species,
                             str_detect(species, 'nucleatum')~species,
                             TRUE~'')) %>% 
  ggplot(aes(x=abs(fc), y=abs(prev.shift), col=colour_group)) + 
    geom_point() + 
    xlab('Absolute generalized fold change') + 
    ylab('Absolute prevalence shift)') + 
    theme_bw() + 
    geom_text_repel(aes(label=highlight), min.segment.length = 0) + 
    scale_colour_manual(values=c('#70737250', '#8BB8E8', 
                                 '#307FE2', '#003DA5'))

# plot CD
ext.measures %>% 
  filter(group=='CD') %>% 
  mutate(species=str_remove(species, 
                            '\\[(ref|meta)_mOTU_v25_[0-9]{5}\\]$')) %>% 
  mutate(fc=abs(fc), prev.shift=abs(prev.shift)) %>% 
  mutate(highlight=case_when(
    fc > 0.75 & prev.shift < 0.2 ~ species,
    prev.shift > 0.5~species,
    fc > 1.4 ~ species,
    TRUE~'')) %>% 
  mutate(unc=str_detect(species, 
                        'Clostridiales species incertae sedis')) %>% 
  ggplot(aes(x=fc, y=prev.shift, col=colour_group)) + 
    geom_point(aes(shape=unc)) + 
    xlab('Absolute generalized fold change') + 
    ylab('Absolute prevalence shift)') + 
    theme_bw() + 
    geom_text_repel(aes(label=highlight), min.segment.length = 0) + 
    scale_colour_manual(values=c('#70737250', '#8BB8E8', 
                                 '#307FE2', '#003DA5')) + 
    scale_shape_manual(values=c(16, 1))

# plot together
g <- sim.measures %>% 
  mutate(fc=abs(fc.sim)) %>% 
  mutate(prev.shift=abs(prev.shift.sim)) %>% 
  transmute(prev.shift, fc, colour_group=case_when(selection~'implanted', 
                                                   TRUE~'background'), 
            group=simulation, type='simulated', 
            selection=as.character(selection)) %>% 
  bind_rows(ext.measures %>% mutate(fc=abs(fc), 
                                    prev.shift=abs(prev.shift),
                                    selection='real')) %>% 
  arrange(p.adj, selection) %>% 
  ggplot(aes(x=fc, y=prev.shift, col=colour_group)) + 
    geom_point() + 
    facet_grid(~group) + 
    xlab('Abs. gFC') + 
    ylab("Abs log2(prevalence shift)") + 
    theme_publication() + 
    scale_colour_manual(values=c('#A8A99E50', '#E4004675', '#70737250', 
                                 '#8BB8E875', '#307FE275', '#003DA575'))
ggsave(g, filename = here('figures', 'figure_reality', 'effect_size_exp.pdf'), 
       useDingbats=FALSE, width = 160, height = 42, units = 'mm')

# show some examples
# get original data
# load(here('files', 'external_data.RData'))
# get sim_all features
temp <- rhdf5::h5read('./simulations/resampling/sim_all.h5', 
                      name='ab4_prev2_rep22')
sim.feat.all <- temp$features
rownames(sim.feat.all) <- temp$feature_names
colnames(sim.feat.all) <- temp$sample_names
sim.feat.all <- prop.table(sim.feat.all, 2)
meta.all <- tibble(Sample_ID=temp$sample_names, 
                   Group=temp$labels, type='sim_all')

temp <- rhdf5::h5read('./simulations/resampling/sim_low.h5',
                      name='ab4_prev2_rep22')
sim.feat.low <- temp$features
rownames(sim.feat.low) <- temp$feature_names
colnames(sim.feat.low) <- temp$sample_names
sim.feat.low <- prop.table(sim.feat.low, 2)
meta.low <- tibble(Sample_ID=temp$sample_names, 
                   Group=-temp$labels, type='sim_low')

# select species
# CRC
sp.1 <- 'Parvimonas micra [ref_mOTU_v25_04287]'
sp.2 <- 'Fusobacterium nucleatum subsp. nucleatum [ref_mOTU_v25_01003]'
# IBD
sp.3 <- 'Clostridiales species incertae sedis [meta_mOTU_v25_12635]'
sp.4 <- '[Ruminococcus] gnavus [ref_mOTU_v25_01594]'#'Proteobacteria sp. [ref_mOTU_v25_00095]'
# all
sp.5 <- 'bact_otu_1823'
sp.6 <- 'bact_otu_11650' #
# low
sp.7 <- 'bact_otu_6585'
sp.8 <- 'bact_otu_3995'

sim.measures %>% 
  mutate(highlight=features %in% c(sp.5, sp.6, sp.7, sp.8)) %>% 
  ggplot(aes(x=abs(fc.sim), y=abs(prev.shift.sim), col=selection)) + 
  geom_point() + 
  xlab('Absolute generalized fold change') + 
  ylab('Absolute prevalence shift)') + 
  theme_bw() +
  facet_grid(~simulation) + 
  scale_colour_manual(values=c('#E4004675', '#FFA30075')) +
  geom_point(data=. %>% filter(highlight), col='black')


df.plot <-
  t(feat.crc[c(sp.1, sp.2),meta.crc$Sample_ID]) %>% 
  as_tibble(rownames='Sample_ID') %>% 
  pivot_longer(-Sample_ID, names_to = 'species', values_to = 'rel.ab') %>% 
  left_join(meta.crc, by='Sample_ID') %>% 
  bind_rows(
    t(feat.cd[c(sp.3, sp.4),meta.cd$Sample_ID]) %>% 
      as_tibble(rownames='Sample_ID') %>% 
      pivot_longer(-Sample_ID, names_to = 'species', values_to = 'rel.ab') %>% 
      left_join(meta.cd, by='Sample_ID')) %>% 
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
                             TRUE~Group)))  %>%
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
  mutate(type=factor(type, levels = c('CRC', 'sim_low', 'CD', 'sim_all'))) %>% 
  ggplot(aes(x=species, y=rel.ab, fill=Group)) + 
    geom_quantile_box(alpha=0.75) + 
    geom_jitter(aes(col=Group), 
                position = position_jitterdodge(jitter.width = 0.2)) + 
    scale_fill_manual(values=c('#70737250', '#307FE2','#307FE2', 
                               '#FFA300', '#734595')) + 
    scale_colour_manual(values=c('#707372', '#307FE2','#307FE2', 
                                 '#FFA300', '#734595')) + 
    theme_publication() + 
    xlab('') + 
    theme(axis.ticks.x = element_blank()) + 
    facet_grid(~type, scales = 'free')
ggsave(g, filename = here('figures', 'figure_reality', 'selected_markers.pdf'),
       useDingbats=FALSE, units = 'mm', height = 70, width = 180)

# all effect sizes for the supplement
var.sum <- var.all %>% 
  filter(simulation!='reference') %>% 
  group_by(simulation, ab, prev, rep, selection) %>% 
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

colour.scheme <- c(
  '#A8C700', '#b4ff16', '#a2ee00', '#A8C700', '#6ca000',
  '#FFCD00', '#ffe476', '#ffdc4e', '#FFCD00', '#c49e00',
  '#FFA300', '#ffd589', '#ffb83b', '#FFA300', '#d88a00',
  '#E40046', '#ff5b8d', '#ff0c57', '#E40046', '#bd003a',
  '#009F4D', '#02ff7d', '#00c660', '#009F4D', '#00783a',
  '#307FE2', '#97bff0', '#649fe9', '#307FE2', '#1d6bce',
  '#8246AF', '#c3a4db', '#9e6cc4', '#8246AF', '#6d3b93'
)

g <- var.sum %>% 
  mutate(group=paste0(ab, '-', prev)) %>% 
  ggplot(aes(x=fc, y=prev.shift, col=group, shape=selection)) + 
    geom_point() + 
    facet_wrap(~simulation) + 
    theme_publication(panel.grid = 'major') + 
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
