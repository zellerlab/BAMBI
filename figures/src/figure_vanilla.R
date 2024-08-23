# ##############################################################################
#
## Figure 2: Results from the vanilla benchmark
#
# ##############################################################################

library("tidyverse")
library("here")
library("ggembl")
library("vroom")

test.colours <- yaml::yaml.load_file(here('files', 'test_colours.yml'))

# ##############################################################################
# Functions

# for (i in c('Zeevi_WGS_sparseDOSSA', 'Zeevi_WGS_SimMSeq', 'Zeevi_WGS_negbin', 
#             'Zeevi_WGS_betabin', 'Zeevi_WGS_all', 'TwinsUK_WGS_all', 
#             'Schirmer_WGS_all', 'HMP_vagina_all', 'HMP_stool_all', 
#             'HMP_skin_all', 'HMP_saliva_all', 'HMP_airways_all')){
#   message(i)
#   .f_read_single_test(fn=paste0("revisions/sim_", i, "_t-test.tsv"),
#                       fn.final=paste0("./figures/figure_vanilla/test_results_combined/res_sim_", i, ".tsv"))
# }

.f_read_single_test <- function(fn, fn.final){
  message(fn)
  res <- vroom::vroom(here('test_results', fn)) %>% 
    separate(group, into=c('group', 'repetition'), 
             sep='_rep') %>% 
    mutate(subset=str_remove(subset, 'subset_')) %>% 
    mutate(subset=as.factor(as.numeric(subset))) %>% 
    mutate(PR = 1-FDR) %>% 
    mutate(FPR=FP/(FP+TN)) %>% 
    mutate(norm=case_when(test=='ZIBSeq-sqrt'~paste0(norm, '-sqrt'), 
                          TRUE~norm)) %>% 
    mutate(norm=case_when(norm=='pass-sqrt'~'sqrt', TRUE~norm)) %>% 
    mutate(test=str_remove(test, '-sqrt$')) %>% 
    select(-problem, -time.running)
  res.combined <- res %>% 
    group_by(group, subset, test, norm, adjust) %>% 
    summarise(precision=mean(PR), sd.precision=sd(PR),
              recall=mean(R), sd.recall=sd(R),
              AUC=mean(auroc), sd.AUC=sd(auroc),
              fdr=mean(FDR), sd.fdr=sd(FDR),
              fpr=mean(FPR), sd.fpr=sd(FPR),
              n.obs=n(),
              .groups='drop') %>% 
    separate(group, into=c('ab','prev'), 
             sep = '_', remove=FALSE) %>% 
    mutate(type=paste0(test, '-', norm))
  # fix ANCOM
  # message(length(unique(res.combined$test)))
  if (unique(res.combined$test)=='ANCOM'){
    res.combined <- res.combined %>% filter(adjust=='none')
  }
  write_tsv(res.combined, file=fn.final, append = file.exists(fn.final))
}

.f_load_and_preprocess <- function(sim){
  
  fn.final <- here('figures', 'figure_vanilla', 'test_results_combined', 
                   paste0('res_', sim, '.tsv'))
  
  fn.res <- list.files(here('test_results'), recursive = TRUE,
                       pattern = paste0(sim, '.+tsv$'))
  if (length(fn.res) == 0 & file.exists(fn.final)){
    df.res <- read_tsv(fn.final, col_types = cols())
  } else if (length(fn.res) == 1 & !file.exists(fn.final)){
    message("+ Loading and combining all results from a single file!")
    df.res <- vroom::vroom(here('test_results', fn.res)) %>% 
      distinct() %>% 
      separate(group, into=c('group', 'repetition'), sep='_rep') %>% 
      mutate(subset=str_remove(subset, 'subset_')) %>% 
      mutate(subset=as.numeric(subset)) %>% 
      mutate(PR = 1-FDR) %>% 
      mutate(FPR=FP/(FP+TN)) %>% 
      mutate(norm=case_when(test=='ZIBSeq-sqrt'~paste0(norm, '-sqrt'), 
                            TRUE~norm)) %>% 
      mutate(test=str_remove(test, '-sqrt$')) %>% 
      select(-problem, -time.running) %>% 
      group_by(group, subset, test, norm, adjust) %>% 
      summarise(precision=mean(PR), sd.precision=sd(PR),
                recall=mean(R), sd.recall=sd(R),
                AUC=mean(auroc), sd.AUC=sd(auroc),
                fdr=mean(FDR), sd.fdr=sd(FDR),
                fpr=mean(FPR), sd.fpr=sd(FPR),
                n.jobs=n(),
                .groups='drop') %>% 
      separate(group, into=c('ab','prev'), 
               sep = '_', remove=FALSE) %>% 
      mutate(type=paste0(test, '-', norm))
    # adjust ANCOM
    df.ancom <- df.res %>% 
      filter(test=='ANCOM') %>% 
      filter(adjust=='none')
    df.res <- df.res %>% 
      filter(test!='ANCOM') %>% 
      bind_rows(df.ancom) %>% 
      bind_rows(df.ancom %>% mutate(adjust='BH')) %>% 
      bind_rows(df.ancom %>% mutate(adjust='BY'))
    write_tsv(df.res, file = fn.final)  
    df.res <- read_tsv(fn.final, col_types=cols())
  } else if (length(fn.res) > 1 & !file.exists(fn.final)){
    message("+ Loading and combining all results from various files!")
    df.res <- map(fn.res, .f = .f_read_single_test, fn.final=fn.final)
    df.res <- read_tsv(fn.final, col_types=cols())
  } else {
    df.res <- read_tsv(fn.final, col_types = cols())
    finished.tests <- unique(df.res$test)
    files.test <- str_remove(fn.res, '.tsv$') %>% 
      str_remove('.*_') %>% 
      intersect(names(test.colours))
    left.over.tests <- setdiff(files.test, finished.tests)
    
    
    if (length(left.over.tests) > 0){
      fn.res.left <- map(left.over.tests, .f = function(x){
        fn.res[str_detect(fn.res, paste0(x, '.tsv'))]}) %>% unlist()
      message("+ Adding results to already existing file for tests: ", 
              paste(left.over.tests, collapse = ', '))
      df.res <- map(fn.res.left, .f = .f_read_single_test, fn.final=fn.final)
      df.res <- read_tsv(fn.final, col_types = cols())
    }
  }
  
  return(df.res)
}

.f_find_norms <- function(df.res){
  norm.rank <- df.res %>% 
    mutate(adjust=case_when(test=='ANCOM'~'BH', TRUE~adjust)) %>% 
    filter(!ab %in% c('ab6', 'ab7')) %>% 
    filter(adjust=='BH') %>% 
    filter(subset >= 100) %>% 
    filter(subset <= 400) %>% 
    select(group, subset, test, norm, recall, fdr, AUC) %>% 
    group_by(test, norm) %>% 
    summarize(n.fdr=mean(fdr > 0.2), 
              m.fdr=mean(fdr),
              m.recall=mean(recall), 
              m.auroc=mean(AUC), .groups='drop') %>% 
    arrange(n.fdr, desc(m.recall)) %>% 
    group_by(test) %>% 
    slice(1) %>% 
    mutate(g=paste0(test, '-', norm))
}

.f_find_test_ranks <- function(df.res, df.norm){
  df.res %>% 
    bind_rows(df.res %>% 
                filter(test=='ANCOM') %>% 
                mutate(adjust='BH')) %>% 
    filter(adjust %in% c('BH', 'none')) %>% 
    mutate(g=paste0(test, '-', norm)) %>% 
    filter(g %in% df.norm$g) %>% 
    filter(subset %in% c(50, 100, 200)) %>%
    filter(!ab %in% c('ab6', 'ab7')) %>% 
    mutate(AUC=case_when(adjust=='none'~AUC, TRUE~NA_real_)) %>% 
    mutate(recall=case_when(adjust=='BH'~recall, TRUE~NA_real_)) %>% 
    mutate(AUC=case_when(!str_detect(group, 'ab1')~AUC, TRUE~NA_real_)) %>% 
    mutate(recall=case_when(!str_detect(group, 'ab1')~recall, TRUE~NA_real_)) %>% 
    mutate(fdr=case_when(adjust=='BH'~fdr, TRUE~NA_real_)) %>% 
    group_by(test) %>% 
    summarise(m.fdr.10=mean(fdr > 0.1, na.rm=TRUE), 
              m.fdr.20=mean(fdr > 0.2, na.rm=TRUE), 
              mean.auroc=mean(AUC, na.rm=TRUE), 
              mean.recall=mean(recall, na.rm=TRUE), 
              mean.fdr=mean(fdr, na.rm=TRUE))#  %>% 
    # mutate(n.fdr=as.numeric(m.fdr > 0.2)) %>%
    # arrange(n.fdr, desc(mean.auroc))
}

norm.levels <- c('rarefy', 'sqrt', 'rarefy-sqrt',
                 'rank',  'clr', 'rclr',
                 'TSS', 'TSS.log', "TSS.arcsin", 
                 'rarefy.TSS',  'rarefy.TSS.log')

# ##############################################################################
# select the right norm for each test, given BH adjustment
df.sim.all <- .f_load_and_preprocess('sim_Zeevi_WGS_all')
norm.rank <- .f_find_norms(df.sim.all)


# compare limma and lm
lm.limma.diff <- df.sim.all %>% 
  filter(test %in% c('lm', 'limma')) %>%
  filter(!str_detect(group, 'ab[67]')) %>% 
  select(group, subset, test, norm, adjust, precision, recall, AUC) %>% 
  pivot_wider(values_from = c(precision, recall, AUC), names_from = test) %>% 
  pivot_longer(cols=c(precision_limma, precision_lm, recall_limma, recall_lm,
                      AUC_limma, AUC_lm)) %>% 
  separate(name, into = c('type', 'test'), sep='_') %>% 
  pivot_wider(values_from = value, names_from = test) %>% 
  mutate(limma=case_when(type=='AUC'~(limma-0.5)/0.5, TRUE~limma)) %>%
  mutate(lm=case_when(type=='AUC'~(lm-0.5)/0.5, TRUE~lm)) %>% 
  filter(!is.na(lm)) %>% 
  group_by(subset, norm, adjust, type) %>% 
  mutate(diff=abs(limma-lm))
lm.limma.diff %>% 
  ggplot(aes(x=as.factor(subset), y=diff, fill=type)) + 
    geom_boxplot() + 
    facet_grid(~adjust)
lm.limma.diff %>% 
  filter(subset>=100) %>% 
  ggplot(aes(x=norm, y=diff, fill=type)) + 
  geom_boxplot() + 
  facet_grid(~adjust)


# compare the different adjustment strategies
df.adjust <- df.sim.all %>% 
  filter(type %in% norm.rank$g) %>% 
  select(test, subset, group, adjust, precision, recall, AUC) %>% 
  mutate(AUC=case_when(AUC < 0.5~0.5, TRUE~AUC)) %>%
  pivot_longer(cols = c(precision, recall, AUC))

df.adjust %>% 
  filter(name=='precision') %>% 
  mutate(fail=value < 0.9) %>% 
  filter(!str_detect(group, 'ab[67]')) %>% 
  group_by(test, adjust, subset) %>% 
  summarise(m=mean(fail)) %>% 
  mutate(adjust=factor(adjust, levels = c('none', 'BH', 'BY'))) %>% 
  mutate(test=factor(test, levels = rev(c('all', names(test.colours))))) %>% 
  mutate(l=sprintf(fmt='%.2f', m)) %>% 
  ggplot(aes(x=as.factor(subset), y=test, fill=m)) + 
    facet_grid(~adjust) +
    geom_tile() + 
    geom_text(aes(label=l), col='white') + 
    xlab('') + ylab('') + 
    theme_publication() +
    scale_fill_gradientn(colors = c('#f3adad', '#8c1515'))

g.adjust <- df.adjust %>% 
  mutate(adjust=factor(adjust, levels = c('none', 'BH', 'BY'))) %>% 
  mutate(test=factor(test, levels = c('all', names(test.colours)))) %>% 
  mutate(name=factor(name, levels=c('precision', 'recall', 'AUC'))) %>% 
  ggplot(aes(x=test, y=value, fill=adjust)) + 
    geom_boxplot() +
    facet_grid(name~., scales = 'free', space = 'free_x') + 
    xlab('') + ylab('') +
    theme_publication(panel.grid = 'major_y') + 
    scale_fill_manual(values=c('none'='#899DA4', 
                               'BH'='#3B9AB2', 
                               'BY'='#F21A00')) + 
    theme(axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5))
ggsave(g.adjust, filename = here('figures', 'figure_vanilla',
                                 'compare_adjustment.pdf'),
       width = 8, height = 4.5, useDingbats=FALSE)

# test ANCOM AUC
ancom.auc <- df.sim.all %>% 
  filter(test==c('ANCOM')) %>% 
  filter(adjust=='none') %>% 
  filter(norm=='pass') %>% 
  filter(!ab %in% c('ab6', 'ab7')) %>% 
  filter(subset < 240) 
wilcox.test(ancom.auc$AUC, mu=0.5, alternative='greater')
t.test(ancom.auc$AUC, mu=0.5, alternative='greater')

# rank tests
ranks.all <- .f_find_test_ranks(df.sim.all, norm.rank)
# line plot
g <- df.sim.all %>% 
  bind_rows(df.sim.all %>% 
              filter(group=='ab4_prev2', test=='ANCOM') %>% 
              mutate(adjust='BH')) %>% 
  filter(group=='ab4_prev2') %>% 
  mutate(g=paste0(test, '-', norm)) %>% 
  filter(g %in% norm.rank$g) %>% 
  select(subset, test, precision, recall, AUC, fpr, adjust) %>% 
  mutate(test=factor(test, levels = names(test.colours))) %>% 
  mutate(fpr=1-fpr) %>% 
  pivot_longer(c(precision, recall, AUC, fpr)) %>% 
  ggplot(aes(x=as.factor(subset), y=value)) + 
  geom_line(aes(group=test, colour=test)) + 
  scale_colour_manual(values=unlist(test.colours)) +
  facet_grid(adjust~name, scales = 'free_y') +
  theme_presentation(panel.grid = 'major_y')
ggsave(g, filename = here('figures', 'figure_vanilla', 'line_plots.pdf'),
       width = 10, height = 7, useDingbats=FALSE)

df.sim.all %>% 
  bind_rows(df.sim.all %>% 
              filter(group=='ab4_prev2', test=='ANCOM') %>% 
              mutate(adjust='BH')) %>% 
  filter(group=='ab4_prev2') %>% 
  mutate(g=paste0(test, '-', norm)) %>% 
  filter(g %in% norm.rank$g) %>% 
  select(subset, test, precision, recall, AUC, fpr, adjust) %>% 
  filter(subset==800, adjust=='none') %>% 
  arrange(AUC)

# compare other norms
norm.data <- df.sim.all %>% 
  bind_rows(df.sim.all %>% 
              filter(test=='ANCOM') %>% 
              mutate(adjust='BH')) %>% 
  filter(adjust%in%c('none', 'BH')) %>% 
  filter(!ab %in% c('ab6', 'ab7')) %>% 
  filter(subset >= 100) %>% 
  filter(subset <= 400) %>% 
  select(group, adjust, subset, norm, test, adjust, precision, recall, AUC) %>% 
  pivot_longer(c(precision, recall, AUC)) %>% 
  mutate(type=paste0(adjust, '-', name)) %>%
  filter(type %in% c('BH-recall', 'BH-precision', 'none-AUC')) %>%
  select(-type) %>% 
  pivot_wider(names_from = norm, values_from = value) %>% 
  pivot_longer(cols=setdiff(norm.levels, c('pass')), 
               names_to = 'norm') %>% 
  filter(!is.na(value)) %>% 
  mutate(change_to_ref=pass-value)

# barchart
g.norm.bar <- norm.data %>% 
  select(group, test, name, norm, adjust, change_to_ref) %>% 
  group_by(test, name, norm, adjust) %>% 
  summarise(m=mean(-change_to_ref, na.rm=TRUE), 
            s1=quantile(-change_to_ref, probs=0.25, na.rm=TRUE),
            s2=quantile(-change_to_ref, probs=0.75, na.rm=TRUE),
            .groups = 'drop') %>% 
  mutate(norm=factor(norm, levels=norm.levels)) %>% 
  mutate(name=factor(name, levels = c('precision', 'recall', 'AUC'))) %>% 
  mutate(test=factor(test, levels=names(test.colours))) %>% 
  ggplot(aes(x=test, y=m, fill=test)) + 
  geom_bar(stat='identity', position = position_dodge()) + 
  facet_grid(name~norm, scale='free_x', space='free') + 
  scale_fill_manual(values=test.colours) +
  theme_publication(panel.grid = 'major_y') + 
  xlab('') + ylab('Change compared to no data transformation') +
  coord_cartesian(ylim=c(-0.3, 0.3)) +
  geom_errorbar(aes(ymin=s1, ymax=s2), width=0.2)
ggsave(g.norm.bar, filename = here('figures', 'figure_vanilla', 
                                   'sim_all_norm.pdf'),
       width = 8, height = 4, useDingbats=FALSE)

# ##############################################################################
# other effect sizes

g.ef <- df.sim.all %>% 
  bind_rows(df.sim.all %>% filter(test=='ANCOM') %>% mutate(adjust='BH')) %>% 
  filter(subset==100) %>% 
  filter(!str_detect(group, 'ab[67]')) %>% 
  mutate(g=paste0(test, '-', norm)) %>% 
  filter(g %in% norm.rank$g) %>% 
  select(ab, adjust, prev, test, precision, recall, AUC, fpr) %>% 
  mutate(fpr=1-fpr) %>% 
  mutate(test=factor(test, levels = names(test.colours))) %>% 
  pivot_longer(c(precision, recall, AUC, fpr)) %>% 
  mutate(type=paste0(adjust, '-', name)) %>%
  filter(type %in% c('BH-recall', 'BH-precision', 'none-AUC', 'none-fpr')) %>%
  select(-type) %>% 
  mutate(name=factor(name, levels=c('precision', 'recall', 'AUC', 'fpr'))) %>% 
  ggplot(aes(x=as.factor(ab), y=value)) + 
    geom_line(aes(group=test, colour=test)) + 
    scale_colour_manual(values=unlist(test.colours)) +
    facet_grid(name~prev) + 
    theme_presentation(panel.grid = 'major_y') + 
    xlab('Abundance scaing') + ylab('')
ggsave(g.ef, filename = here('figures', 'figure_vanilla', 'effect_size.pdf'),
       width = 8, height = 5, useDingbats=FALSE)

# ##############################################################################
# other simulations
# df.sim.MMH <- .f_load_and_preprocess('sim_Zeevi_WGS_MMH')
df.sim.Weiss <- .f_load_and_preprocess('sim_Zeevi_WGS_Weiss')
# df.sim.dir <- .f_load_and_preprocess('sim_Zeevi_WGS_dir')
df.sim.beta <- .f_load_and_preprocess('sim_Zeevi_WGS_betabin')
df.sim.neg <- .f_load_and_preprocess('sim_Zeevi_WGS_negbin')
# df.sim.nbc <- .f_load_and_preprocess('sim_Zeevi_WGS_nbc')
df.sim.spd <- .f_load_and_preprocess('sim_Zeevi_WGS_sparseDOSSA')
df.sim.sim <- .f_load_and_preprocess('sim_Zeevi_WGS_SimMSeq')

# norm.rank.MMH <- .f_find_norms(df.sim.MMH)
norm.rank.Weiss <- .f_find_norms(df.sim.Weiss)
# norm.rank.dir <- .f_find_norms(df.sim.dir)
norm.rank.beta <- .f_find_norms(df.sim.beta)
norm.rank.neg <- .f_find_norms(df.sim.neg)
# norm.rank.nbc <- .f_find_norms(df.sim.nbc)
norm.rank.spd <- .f_find_norms(df.sim.spd)
norm.rank.sim <- .f_find_norms(df.sim.sim)

map(list(#norm.rank.MMH, 
  norm.rank.Weiss, # norm.rank.dir, 
  norm.rank.beta,
  norm.rank.neg, #norm.rank.nbc, 
  norm.rank.spd,
  norm.rank.sim), nrow)

# ranks.MMH <- .f_find_test_ranks(df.sim.MMH, norm.rank.MMH) %>% 
#   mutate(type='MMH')
ranks.Weiss <- .f_find_test_ranks(df.sim.Weiss, norm.rank.Weiss) %>% 
  mutate(type='Weiss')
# ranks.dir <- .f_find_test_ranks(df.sim.dir, norm.rank.dir) %>% 
#   mutate(type='dir')
ranks.beta <- .f_find_test_ranks(df.sim.beta, norm.rank.beta) %>% 
  mutate(type='betabin')
ranks.neg <- .f_find_test_ranks(df.sim.neg, norm.rank.neg) %>% 
  mutate(type='negbin')
# ranks.nbc <- .f_find_test_ranks(df.sim.nbc, norm.rank.nbc) %>% 
#   mutate(type='nbc')
ranks.spd <- .f_find_test_ranks(df.sim.spd, norm.rank.spd) %>% 
  mutate(type='sparseDOSSA')
ranks.sim <- .f_find_test_ranks(df.sim.sim, norm.rank.sim) %>% 
  mutate(type='SimMSeq')
ranks.all <- .f_find_test_ranks(df.sim.all %>% filter(test!='LDM'), norm.rank) %>% 
  mutate(type='Zeevi_WGS')

# ##############################################################################
# other datasets comparison
df.TwinsUK <- .f_load_and_preprocess('sim_TwinsUK_WGS_all')
df.Schirmer <- .f_load_and_preprocess('sim_Schirmer_WGS_all')
df.airways <- .f_load_and_preprocess('sim_HMP_airways')
df.saliva <- .f_load_and_preprocess('sim_HMP_saliva')
df.skin <- .f_load_and_preprocess('sim_HMP_skin')
df.vagina <- .f_load_and_preprocess('sim_HMP_vagina')
df.stool <- .f_load_and_preprocess('sim_HMP_stool')

norm.rank.twins <- .f_find_norms(df.TwinsUK)
norm.rank.schirmer <- .f_find_norms(df.Schirmer)
norm.rank.airways <- .f_find_norms(df.airways)
norm.rank.saliva <- .f_find_norms(df.saliva)
norm.rank.skin <- .f_find_norms(df.skin)
norm.rank.stool <- .f_find_norms(df.stool)
norm.rank.vagina <- .f_find_norms(df.vagina %>% 
                                    mutate(subset=case_when(subset==94~100, 
                                                            TRUE~subset)))

map(list(norm.rank.twins, norm.rank.schirmer, norm.rank.airways, 
         norm.rank.stool, norm.rank.saliva, norm.rank.skin, norm.rank.vagina),nrow)

rank.tests.twins <- .f_find_test_ranks(df.TwinsUK, norm.rank)
rank.tests.schirmer <- .f_find_test_ranks(df.Schirmer, norm.rank)
rank.tests.airways <- .f_find_test_ranks(df.airways, norm.rank)
rank.tests.saliva <- .f_find_test_ranks(df.saliva, norm.rank)
rank.tests.stool <- .f_find_test_ranks(df.stool, norm.rank)
rank.tests.skin <- .f_find_test_ranks(df.skin, norm.rank)
rank.tests.vagina <- .f_find_test_ranks(df.vagina, norm.rank)



df.heat <- bind_rows(#ranks.MMH, 
  ranks.Weiss, 
  ranks.spd, 
  ranks.beta, #ranks.dir, 
  ranks.neg, #ranks.nbc, 
  ranks.sim,
                    ranks.all %>% mutate(type='Zeevi'), 
                    rank.tests.twins %>% mutate(type='TwinsUK'), 
                    rank.tests.schirmer %>% mutate(type='Schirmer'), 
                    rank.tests.airways %>% mutate(type='HMP_airways'),  
                    rank.tests.saliva %>% mutate(type='HMP_saliva'),  
                    rank.tests.skin %>% mutate(type='HMP_skin'), 
                    rank.tests.stool %>% mutate(type='HMP_stool'), 
                    rank.tests.vagina %>% mutate(type='HMP_vagina')) %>% 
  mutate(n.fdr=as.numeric(m.fdr.10 > 0.1)) %>% 
  filter(test!='LDM') %>% 
  group_by(type) %>% 
  mutate(rank.auroc=case_when(n.fdr==1~NA_real_, TRUE~mean.auroc)) %>% 
  mutate(test.rank=rank(desc(rank.auroc))) %>% 
  ungroup() %>% 
  mutate(test.rank=case_when(n.fdr==1~NA_real_, TRUE~test.rank)) %>% 
  mutate(test=factor(test, levels = ranks.all %>% 
                       mutate(ranking=mean.auroc) %>% 
                       arrange(m.fdr.10, desc(ranking)) %>% 
                       pull(test) %>% rev)) %>% 
  mutate(type=factor(
    type, levels = c('MMH', 'Weiss', 'betabin', 'dir', 'negbin', 'nbc', 
                     'SimMSeq', 'sparseDOSSA',
                     'Zeevi', 'Schirmer', 'TwinsUK', 'HMP_stool', 
                     'HMP_skin', 'HMP_saliva', 'HMP_airways', 
                     'HMP_vagina')))
g.heat <- df.heat %>% 
  mutate(mean.auroc=case_when(n.fdr==0~mean.auroc, n.fdr==1~0.4)) %>% 
  mutate(cell.label=
           case_when(test.rank < 4~paste0(test.rank, ' - ', 
                                          sprintf(fmt='%.2f', mean.auroc)),
                     is.na(test.rank)~'',
                     TRUE~sprintf(fmt='%.2f', mean.auroc))) %>% 
  ggplot(aes(x=type, y=test, fill=mean.auroc)) + 
  geom_tile() + 
  scale_fill_gradientn(colours=c('#D0DEBB', '#6CC24A', '#18974C', 
                                 '#007B53', '#0A5032'),
                       limits=c(0.5, 1), na.value = '#D0D0CE') +
  theme_publication() + 
  geom_text(aes(label=cell.label)) + 
  xlab('') + ylab('')

g.heat.bw <- df.heat %>% 
  mutate(cell.label=
           case_when(test.rank < 4~'',
                     is.na(test.rank)~sprintf(fmt='%.2f', mean.auroc),
                     TRUE~'')) %>% 
  ggplot(aes(x=type, y=test, fill=mean.auroc)) + 
    geom_tile() + 
    scale_fill_gradientn(colours=c('#D0D0CE', '#A8A99E', '#707372', 
                                   '#54585A', '#373A36'),
                         limits=c(0.5, 1), na.value = 'white') +
    theme_publication() + 
    geom_text(aes(label=cell.label)) + 
    xlab('') + ylab('')


ggsave(g.heat, filename = here('figures', 'figure_vanilla', 
                               'comparison_everything.pdf'),
       width = 12, height = 6, useDingbats=FALSE)
ggsave(g.heat.bw, filename = here('figures', 'figure_vanilla', 
                                  'comparison_everything_bw.pdf'),
       width = 12, height = 6, useDingbats=FALSE)


# correlation
g <- df.sim.all %>% 
  mutate(adjust=case_when(test=='ANCOM'~'BH', TRUE~adjust)) %>% 
  filter(adjust=='BH') %>% 
  filter(type %in% norm.rank$g) %>% 
  # filter(prev=='prev2') %>% 
  # mutate(group=str_remove(group, '_prev[0-9]')) %>% 
  select(group, subset, test, precision, recall, AUC) %>% 
  pivot_longer(c(precision, recall, AUC)) %>% 
  full_join(df.TwinsUK %>% 
              mutate(adjust=case_when(test=='ANCOM'~'BH', TRUE~adjust)) %>% 
              filter(adjust=='BH') %>% 
              filter(type %in% norm.rank$g) %>% 
              select(group, subset, test, precision, recall, AUC) %>% 
              pivot_longer(c(precision, recall, AUC), 
                           values_to='TwinsUK'),
            by=c("group", "subset", "test", "name")) %>% 
  full_join(df.Schirmer %>% 
              mutate(adjust=case_when(test=='ANCOM'~'BH', TRUE~adjust)) %>% 
              filter(adjust=='BH') %>% 
              filter(type %in% norm.rank$g) %>% 
              select(group, subset, test, precision, recall, AUC) %>% 
              pivot_longer(c(precision, recall, AUC), 
                           values_to='Schirmer'),
            by=c("group", "subset", "test", "name"))  %>% 
  pivot_longer(c(TwinsUK, Schirmer), names_to = 'dataset', 
               values_to = 'value_dataset') %>% 
  mutate(test=factor(test, levels = names(test.colours))) %>% 
  mutate(value=case_when(name=='AUC'~(value-0.5)/0.5, TRUE~value)) %>% 
  mutate(value_dataset=case_when(name=='AUC'~(value_dataset-0.5)/0.5, 
                                 TRUE~value_dataset)) %>% 
  mutate(value=abs(value), value_dataset=abs(value_dataset)) %>%
  mutate(name=factor(name, levels = c('precision', 'recall', 'AUC'))) %>% 
  ggplot(aes(x=value, y=value_dataset, col=test)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_point() + 
  scale_colour_manual(values=unlist(test.colours)) + 
  facet_grid(dataset~name) + 
  theme_publication(panel.grid = 'major') +
  xlab('Value from Zeevi WGS') + 
  ylab('Value from other dataset')
ggsave(g, filename = here('figures', 'figure_vanilla', 'correlation_WGS.pdf'),
       width = 8, height = 5, useDingbats=FALSE)


# ##############################################################################
# line plots for other datasets
.f_get_interesting_part <- function(df.res, df.norm, d){
  if (all(is.na(df.res$prev))){
    g <- 'ab4'
  } else {
    g <- 'ab4_prev2'
  }
  df.res %>% 
    bind_rows(df.res %>% 
                filter(group==g, test=='ANCOM') %>% 
                mutate(adjust='BH')) %>%
    filter(group==g) %>% 
    mutate(dataset=d) %>%
    mutate(g=paste0(test, '-', norm)) %>% 
    filter(g %in% df.norm$g) %>% 
    select(subset, test, precision, recall, AUC, fpr, adjust, dataset) %>% 
    mutate(test=factor(test, levels = names(test.colours))) %>% 
    mutate(fpr=1-fpr) %>% 
    pivot_longer(c(precision, recall, AUC, fpr)) %>% 
    mutate(type=paste0(name, '-', adjust)) %>% 
    filter(type %in% c('precision-BH', 'recall-BH', 
                       'AUC-none', 'fpr-none')) %>% 
    select(-type) 
}

df.plot <- .f_get_interesting_part(df.sim.all, norm.rank, 'Zeevi_WGS') %>% 
  bind_rows(.f_get_interesting_part(df.Schirmer, 
                                    norm.rank.schirmer, 'Schirmer_WGS')) %>% 
  bind_rows(.f_get_interesting_part(df.TwinsUK, 
                                    norm.rank.twins, 'TwinsUK_WGS')) %>% 
  # bind_rows(.f_get_interesting_part(df.sim.MMH, 
                                    # norm.rank.MMH, 'MMH')) %>% 
  bind_rows(.f_get_interesting_part(df.sim.Weiss,
                                    norm.rank.Weiss, 'Weiss')) %>% 
  bind_rows(.f_get_interesting_part(df.sim.sim,
                                    norm.rank.sim, 'SimMSeq')) %>%
  bind_rows(.f_get_interesting_part(df.sim.beta, 
                                    norm.rank.beta, 'betabin')) %>% 
  bind_rows(.f_get_interesting_part(df.sim.neg, 
                                    norm.rank.neg, 'negbin')) %>% 
  # bind_rows(.f_get_interesting_part(df.sim.nbc, 
                                    # norm.rank.nbc, 'nbc')) %>% 
  bind_rows(.f_get_interesting_part(df.sim.spd, 
                                    norm.rank.spd, 'sparseDOSSA')) %>% 
  bind_rows(.f_get_interesting_part(df.stool, 
                                    norm.rank.stool, 'HMP_stool')) %>% 
  bind_rows(.f_get_interesting_part(df.skin, 
                                    norm.rank.skin, 'HPM_skin')) %>% 
  bind_rows(.f_get_interesting_part(df.saliva, 
                                    norm.rank.saliva, 'HMP_saliva')) %>% 
  bind_rows(.f_get_interesting_part(df.vagina, 
                                    norm.rank.vagina, 'HMP_vagina')) %>% 
  bind_rows(.f_get_interesting_part(df.airways, 
                                    norm.rank.airways, 'HMP_airways')) 

table(df.plot$subset)
df.plot <- df.plot %>% 
  mutate(subset=case_when(subset==94~100, subset==178~100, subset==184~200,
                          subset==188~200, TRUE~subset))

g.lines.all <- df.plot %>% 
  ggplot(aes(x=as.factor(subset), y=value)) + 
  geom_line(aes(group=test, colour=test)) + 
  scale_colour_manual(values=unlist(test.colours)) +
  facet_grid(dataset~name) + 
  theme_presentation(panel.grid = 'major_y')
ggsave(g.lines.all, filename = here('figures', 'figure_vanilla', 
                                    'line_plots_all.pdf'),
       width = 16, height = 25, useDingbats=FALSE)

# similarity/correlation between Zeevi and other flavors overall?








# ##############################################################################
# func/tax comparison
df.kegg <- .f_load_and_preprocess('sim_Zeevi_KEGG_all')
norm.kegg <- .f_find_norms(df.kegg)


df.comp <- df.kegg %>% 
  filter(adjust=='BH') %>% 
  filter(type %in% norm.kegg$g) %>% 
  mutate(type='func') %>% 
  bind_rows(df.sim.all %>% 
              filter(adjust=='BH') %>% 
              filter(type %in% norm.rank$g) %>% 
              mutate(type='tax')) %>% 
  select(group, subset, test, precision, recall, AUC, type) %>% 
  pivot_longer(c(precision, recall, AUC)) %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  filter(!is.na(func)) %>% 
  pivot_longer(cols = c(tax, func), names_to = 'type') 
  
g.comp.tax.func <- df.comp %>%
  mutate(type=factor(type, levels=c('tax', 'func'))) %>% 
  mutate(type2='all') %>% 
  mutate(test='all') %>% 
  bind_rows(df.comp %>% 
              mutate(type=factor(type, levels=c('tax', 'func'))) %>% 
              mutate(type2='single')) %>% 
  mutate(test=factor(test, levels = c('all', names(test.colours)))) %>% 
  mutate(name=factor(name, levels = c('precision', 'recall', 'AUC'))) %>% 
  ggplot(aes(x=test, y=value, fill=type)) + 
    geom_boxplot() +
    facet_grid(name~type2, scales='free', space='free_x') + 
    scale_fill_manual(values=c('tax'='#009B76', 'func'='#007C92')) + 
    theme_publication(panel.grid = 'major_y') + 
    xlab('') + ylab('')
ggsave(g.comp.tax.func, filename = './figures/misc/compare_tax_func_tests.pdf',
       width = 6.5, height = 4.5, useDingbats=FALSE)








# ##############################################################################

# overall runtime?
sims <- c('sim_HMP_airways', 'sim_HMP_saliva', 'sim_HMP_skin', 'sim_HMP_vagina',
          'sim_Schirmer_WGS', 'sim_TwinsUK_WGS', 'sim_Zeevi_KEGG_all', 
          'sim_Zeevi_WGS_all', 'sim_Zeevi_WGS_others')
time.total <- list()
for (sim in sims){
  message(sim)
  fn.res <- list.files(here('test_results'), recursive = TRUE,
                       pattern = paste0(sim, '.+tsv$'))
  time.sim <- map(fn.res, .f=function(x){
    read_tsv(here('test_results', sim, x), col_types = cols()) %>% 
      select(test, job.id, time.running) %>% 
      distinct()}) %>% 
    bind_rows() %>%
    group_by(test) %>% 
    summarise(time=sum(time.running)) %>% 
    mutate(simulation=sim)
  time.total[[sim]] <- time.sim
}

time.total <- bind_rows(time.total) %>% 
  mutate(time_hours=time/60) %>% mutate(cost=time_hours*0.06) %>% 
  mutate(test=case_when(test=='ZIBSeq-sqrt'~'ZIBSeq', TRUE~test)) %>% 
  group_by(test, simulation) %>% 
  summarise(time=sum(time), time_hours=sum(time_hours), cost=sum(cost))

order <- time.total %>% 
  filter(simulation=='sim_Zeevi_WGS_all') %>% 
  arrange(cost)

g1 <- time.total %>% 
  filter(simulation=='sim_Zeevi_WGS_all') %>% 
  mutate(test=factor(test, levels = order$test)) %>% 
  ggplot(aes(x=test, y=cost, fill=test)) + 
    geom_bar(stat='identity') + 
    xlab('') + ylab('Compute cost on SCG [$]') + 
    theme_presentation() + 
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
    scale_fill_manual(values = test.colours, guide=FALSE)
ggsave(g1, filename = './figures/misc/cost_Zeevi_WGS_all.pdf', 
       width = 6, height = 4, useDingbats=FALSE)
all <- time.total %>% 
  group_by(test) %>% 
  summarise(cost=sum(cost))
g2 <- all %>% 
  mutate(test=factor(test, levels = order$test)) %>% 
  ggplot(aes(x=test, y=cost, fill=test)) + 
    geom_bar(stat='identity') + 
    xlab('') + ylab('Compute cost on SCG [$]') + 
    theme_presentation() + 
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
    scale_fill_manual(values = test.colours, guide=FALSE) + 
    geom_text(aes(label=sprintf('%.f$', cost)), angle=90, hjust=-0.3) + 
    ylim(0, 3000000)
ggsave(g2, filename = './figures/misc/cost_all.pdf', 
       width = 6, height = 4, useDingbats=FALSE)


#












# use the BH correction only for main figures





norm.data %>% 
  group_by(group, subset, test, norm, adjust) %>% 
  summarize(overall=sum(change_to_ref), .groups = 'drop') %>% 
  group_by(test, norm, adjust) %>%
  summarise(m=mean(overall), .groups='drop') %>% 
  arrange(m) %>% 
  group_by(test, adjust) %>% 
  slice(1) %>% View



g <- df.sim.all %>% 
  filter(group=='ab4_prev2') %>% 
  mutate(g=paste0(test, '-', norm, '-', adjust)) %>% 
  filter(g %in% norm.rank$g) %>% 
  select(subset, test, precision, adjust, recall, AUC, fdr) %>% 
  mutate(test=factor(test, levels = names(test.colours))) %>% 
  pivot_longer(c(precision, recall, AUC, fdr)) %>% 
  mutate(value=case_when(name!='fdr'~value, TRUE~log10(value+1e-03))) %>% 
  ggplot(aes(x=as.factor(subset), y=value)) + 
    geom_line(aes(group=test, colour=test)) + 
    scale_colour_manual(values=unlist(test.colours)) +
    facet_wrap(name~adjust, scales = 'free_y') + 
    theme_presentation(panel.grid = 'major_y')

        
ranks.test <- df.sim.all %>% 
  filter(!ab %in% c('ab6', 'ab7')) %>% 
  filter(adjust=='BH') %>% 
  mutate(g=paste0(test, '-', norm)) %>% 
  filter(g %in% norm.rank$g) %>% 
  filter(subset > 30 & subset < 300) %>%
  group_by(type) %>% 
  summarise(m.fdr=mean(fdr), m.auc=mean(AUC), n.fdr=sum(fdr > 0.1), 
            .groups = 'drop') %>% 
  arrange(n.fdr, desc(m.auc)) %>%
  # arrange(m.fdr) %>% 
  separate(type, into=c('test', 'norm'), sep='-')

g.line <- df.sim.all %>% 
  filter(group=='ab4_prev1') %>% 
  filter(adjust=='BH') %>% 
  mutate(g=paste0(test, '-', norm)) %>% 
  filter(g %in% norm.rank$g) %>% 
  mutate(test=factor(test, levels = ranks.test$test)) %>% 
  select(test, subset, fdr, sd.fdr, AUC, sd.AUC, recall, sd.recall) %>% 
  mutate(AUC=case_when(AUC < 0.5~1-AUC, TRUE~AUC)) %>% 
  mutate(AUC=(AUC-0.5)/0.5, sd.AUC=(sd.AUC/0.5)) %>% 
  pivot_longer(cols=c(fdr, AUC, recall)) %>%
  mutate(value.up=case_when(name=='fdr'~value+sd.fdr, 
                            name=='recall'~value+sd.recall,
                            name=='AUC'~value+sd.AUC)) %>% 
  mutate(value.down=case_when(name=='fdr'~value-sd.fdr, 
                              name=='recall'~value-sd.recall,
                              name=='AUC'~value-sd.AUC)) %>% 
  ggplot(aes(x=as.factor(subset), y=value)) + 
    geom_line(aes(group=name, colour=name)) + 
    geom_ribbon(aes(group=name, ymin=value.down, 
                    ymax=value.up, fill=name), alpha=0.25) +
    facet_grid(.~test) + 
    coord_cartesian(ylim = c(0, 1)) + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1)) + 
    xlab("Sample Size") + ylab("False discovery rate") + 
    scale_colour_manual(values=c('#307FE2', '#E40046', '#18974C')) + 
    scale_fill_manual(values=c('#307FE2', '#E40046', '#18974C')) + 
    geom_hline(yintercept = 0.1, colour='black', linetype=3)
  

fdr.log <- df.sim.all %>% 
  filter(group=='ab4_prev2') %>% 
  filter(adjust=='BH') %>% 
  mutate(g=paste0(test, '-', norm)) %>% 
  filter(g %in% norm.rank$g) %>% 
  select(subset, test, precision, recall, AUC, fdr) %>% 
  mutate(test=factor(test, levels = names(test.colours))) %>% 
  pivot_longer(c(precision, recall, AUC, fdr)) %>% 
  filter(name=='fdr') %>% 
  mutate(value=value+1e-03) %>% 
  ggplot(aes(x=as.factor(subset), y=value)) + 
    geom_line(aes(group=test, colour=test)) +
    scale_colour_manual(values=unlist(test.colours)) +
    theme_publication(panel.grid='major_y') +
    scale_y_log10() + annotation_logticks(sides='l')

ggsave(g, filename = here('figures', 'figure_vanilla', 'sim_all_combined.pdf'),
       width = 6, height = 5, useDingbats=FALSE)
ggsave(g.line, 
       filename = here('figures', 'figure_vanilla', 'sim_all_lines.pdf'),
       width = 12, height = 5, useDingbats=FALSE)
ggsave(fdr.log,
       filename = here('figures', 'figure_vanilla', 'sim_all_fdr_log10.pdf'),
       width = 6, height = 5, useDingbats=FALSE)


# how does it look with increasing effect sizes?
g.effect <- df.sim.all %>% 
  filter(subset==200) %>% 
  mutate(g=paste0(test, '-', norm)) %>% 
  filter(g %in% norm.rank$g) %>% 
  select(ab, prev, test, precision, recall, AUC, fdr) %>% 
  mutate(test=factor(test, levels = names(test.colours))) %>% 
  pivot_longer(c(precision, recall, AUC, fdr)) %>% 
  mutate(value=case_when(name!='fdr'~value, TRUE~log10(value+1e-03))) %>% 
  mutate(lt=test=='metagenomeSeq2') %>% 
  ggplot(aes(x=as.factor(ab), y=value)) + 
  geom_line(aes(group=test, colour=test, lty=lt)) + 
  scale_colour_manual(values=unlist(test.colours)) +
  facet_wrap(prev~name, scales = 'free_y') + 
  theme_presentation(panel.grid = 'major_y')
ggsave(g.effect, filename = here('figures', 'figure_vanilla', 
                                 'sim_all_effect_size.pdf'),
       width = 10, height = 10, useDingbats=FALSE)

# ##############################################################################
# Comparison to other simualtions

# other 



g.heat <- bind_rows(ranks.MMH, ranks.Weiss,ranks.spd, 
          ranks.beta,  ranks.dir, ranks.neg,
          ranks.all) %>% 
  mutate(mean.auroc=case_when(n.fdr==0~mean.auroc, n.fdr==1~0.5)) %>% 
  group_by(type) %>% 
  mutate(test.rank=rank(desc(mean.auroc))) %>% 
  ungroup() %>% 
  mutate(test.rank=case_when(n.fdr==1~NA_real_, TRUE~test.rank)) %>% 
  mutate(test=factor(test, levels = rev(names(test.colours)))) %>% 
  mutate(type=factor(type, levels = c('MMH', 'Weiss', 'sparseDOSSA', 
                                      'betabin', 'dir', 'negbin', 
                                      'sim_all'))) %>% 
  mutate(cell.label=
           case_when(test.rank < 4~paste0(test.rank, ' - ', 
                                          sprintf(fmt='%.2f', mean.auroc)),
                              is.na(test.rank)~'',
                              TRUE~sprintf(fmt='%.2f', mean.auroc))) %>% 
  ggplot(aes(x=type, y=test, fill=mean.auroc)) + 
    geom_tile() + 
    scale_fill_gradientn(colours=c("white", '#D0D0CE', '#A8A99E', 
                                   '#707372', '#54585A', '#373A36'), 
                         limits=c(0.5, 1)) +
    theme_publication() + 
    geom_text(aes(label=cell.label)) + 
    xlab('') + ylab('')
ggsave(g.heat, filename = here('figures', 'figure_vanilla', 
                               'comparison_ranks.pdf'),
       width = 7, height = 7, useDingbats=FALSE)

# line plot for other methods as well
g <- bind_rows(
  df.sim.MMH %>% filter(group=='ab4') %>%  
    mutate(g=paste0(test, '-', norm)) %>% 
    filter(g %in% norm.rank.MMH$g) %>% 
    mutate(type='MMH'),
  df.sim.Weiss %>% filter(group=='ab4') %>%  
    mutate(g=paste0(test, '-', norm)) %>% 
    filter(g %in% norm.rank.Weiss$g) %>% 
    mutate(type='Weiss'),
  df.sim.dir %>% filter(group=='ab4') %>%  
    mutate(g=paste0(test, '-', norm)) %>% 
    filter(g %in% norm.rank.dir$g) %>% 
    mutate(type='dirichlet'),
  df.sim.neg %>% filter(group=='ab4') %>%  
    mutate(g=paste0(test, '-', norm)) %>% 
    filter(g %in% norm.rank.neg$g) %>% 
    mutate(type='negbin'),
  df.sim.spd %>% filter(group=='ab4') %>%  
    mutate(g=paste0(test, '-', norm)) %>% 
    filter(g %in% norm.rank.spd$g) %>% 
    mutate(type='sparseDOSSA'),
  df.sim.all %>% filter(group=='ab4_prev2') %>%  
    mutate(g=paste0(test, '-', norm)) %>% 
    filter(g %in% norm.rank$g) %>% 
    mutate(type='implantation'),
  df.sim.beta %>% filter(group=='ab4') %>%
    mutate(g=paste0(test, '-', norm)) %>% 
    filter(g %in% norm.rank.beta$g) %>% 
    mutate(type='betabin')) %>% 
  select(subset, test, precision, recall, AUC, type) %>% 
  mutate(test=factor(test, levels = names(test.colours))) %>% 
  pivot_longer(c(precision, recall, AUC)) %>% 
  mutate(name=factor(name, levels = c('precision', 'recall', 'AUC'))) %>% 
  mutate(lt=test=='metagenomeSeq2') %>% 
  # mutate(value=case_when(name!='fdr'~value, TRUE~log10(value+1e-03))) %>% 
  ggplot(aes(x=as.factor(subset), y=value)) + 
    geom_line(aes(group=test, colour=test, linetype=lt)) + 
    scale_colour_manual(values=unlist(test.colours)) +
    facet_grid(type~name) + 
    theme_presentation(panel.grid = 'major_y')
ggsave(g, filename = here('figures', 'figure_vanilla', 'comparison_lins.pdf'),
       width = 12, height = 12, useDingbats=FALSE)
