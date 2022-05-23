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

.f_load_and_preprocess <- function(sim){
    
    fn.final <- here('figures', 'figure_vanilla', paste0('res_', sim, '.tsv'))
    
    if (file.exists(fn.final)){
        df.res <- read_tsv(fn.final)
    } else {
        fn.res <- list.files(here('test_results'), recursive = TRUE,
                             pattern = paste0(sim, '.+tsv$'))
        if (length(fn.res)==0){
            stop("Results have not yet been retrieved!")
        } else if (length(fn.res)==1) {
            df.res <- vroom::vroom(here('test_results', fn.res)) %>% 
                separate(group, into=c('group', 'repetition'), sep='_rep') %>% 
                mutate(subset=str_remove(subset, 'subset_')) %>% 
                mutate(subset=as.factor(as.numeric(subset))) %>% 
                mutate(PR = 1-FDR) %>% 
                mutate(norm=case_when(test=='ZIBSeq-sqrt'~'sqrt', 
                                      TRUE~norm)) %>% 
                mutate(test=str_remove(test, '-sqrt$')) %>% 
                select(-job.id, -problem, -time.running) %>% 
                group_by(group, subset, test, norm) %>% 
                summarise(precision=mean(PR), sd.precision=sd(PR),
                          recall=mean(R), sd.recall=sd(R),
                          AUC=mean(auroc), sd.AUC=sd(auroc),
                          fdr=mean(FDR), sd.fdr=sd(FDR),
                          .groups='drop') %>% 
                separate(group, into=c('ab','prev'), 
                         sep = '_', remove=FALSE) %>% 
                mutate(type=paste0(test, '-', norm))
            write_tsv(df.res, file = fn.final)
        } else {
            if (any(str_detect(fn.res, 'fisher'))){
                fn.res <- fn.res[-which(str_detect(fn.res, 'fisher'))]
            }
            df.res <- map(fn.res, .f = function(x){
                res <- vroom::vroom(here('test_results', x)) %>% 
                    separate(group, into=c('group', 'repetition'), 
                             sep='_rep') %>% 
                    mutate(subset=str_remove(subset, 'subset_')) %>% 
                    mutate(subset=as.factor(as.numeric(subset))) %>% 
                    mutate(PR = 1-FDR) %>% 
                    mutate(norm=case_when(test=='ZIBSeq-sqrt'~'sqrt', 
                                          TRUE~norm)) %>% 
                    mutate(test=str_remove(test, '-sqrt$')) %>% 
                    select(-job.id, -problem, -time.running) %>% 
                    group_by(group, subset, test, norm) %>% 
                    summarise(precision=mean(PR), sd.precision=sd(PR),
                              recall=mean(R), sd.recall=sd(R),
                              AUC=mean(auroc), sd.AUC=sd(auroc),
                              fdr=mean(FDR), sd.fdr=sd(FDR),
                              .groups='drop') %>% 
                    separate(group, into=c('ab','prev'), 
                             sep = '_', remove=FALSE) %>% 
                    mutate(type=paste0(test, '-', norm))
                write_tsv(res, file=fn.final, append = file.exists(fn.final))
            })
        }
    }
    return(df.res)
}

.f_find_norms <- function(df.res){
  norm.rank <- df.res %>% 
    filter(!ab %in% c('ab6', 'ab7')) %>% 
    filter(subset >= 100) %>% 
    filter(subset <= 400) %>% 
    select(group, subset, test, norm, recall, fdr, AUC) %>% 
    group_by(test, norm) %>% 
    summarize(n.fdr=mean(fdr > 0.2), 
              m.recall=mean(recall), 
              m.auroc=mean(AUC)) %>% 
    arrange(n.fdr, desc(m.recall)) %>% 
    group_by(test) %>% 
    slice(1) %>% 
    mutate(g=paste0(test, '-', norm))
}

df.sim.all <- .f_load_and_preprocess('sim_all')

norm.colours <- c('pass'='darkgrey',
                  'rank'='darkgrey',
                  'clr'='#A1D99B', 'rclr'='#006D2C',
                  'TSS'='#C6DBEF', 'TSS.log'='#6BAED6', 
                  "TSS.arcsin"='#08519C',
                  'rarefy'='#FB6A4A', 
                  'rarefy.TSS'='#BCBDDC', 
                  'rarefy.TSS.log'='#807DBA',
                  'sqrt'='#A0522D')

norm.data <- df.sim.all %>% 
  filter(!ab %in% c('ab6', 'ab7')) %>% 
  filter(subset >= 100) %>% 
  filter(subset <= 400) %>% 
  select(group, subset, norm, test, precision, recall, AUC) %>% 
  pivot_longer(c(precision, recall, AUC)) %>% 
  pivot_wider(names_from = norm, values_from = value) %>% 
  pivot_longer(cols=setdiff(names(norm.colours), c('pass')), 
               names_to = 'norm') %>% 
  filter(!is.na(value)) %>% 
  mutate(change_to_ref=pass-value)

# barchart
g.norm.bar <- norm.data %>% 
  select(group, test, name, norm, change_to_ref) %>% 
  group_by(test, name, norm) %>% 
  summarise(m=mean(-change_to_ref), s1=quantile(-change_to_ref, probs=0.25),
            s2=quantile(-change_to_ref, probs=0.75),
            .groups = 'drop') %>% 
  mutate(norm=factor(norm, levels=names(norm.colours))) %>% 
  mutate(name=factor(name, levels = c('precision', 'recall', 'AUC'))) %>% 
  mutate(test=factor(test, levels=names(test.colours))) %>% 
  ggplot(aes(x=test, y=m, fill=test)) + 
    geom_bar(stat='identity') + 
    facet_grid(name~norm, scale='free_x', space='free_x') + 
    scale_fill_manual(values=test.colours) + 
    theme_publication(panel.grid = 'major_y') + 
    xlab('') + ylab('Change compared to no normalization') +
    coord_cartesian(ylim=c(-0.3, 0.3)) +
    geom_errorbar(aes(ymin=s1, ymax=s2), width=0.2)

# heatmap
g.norm.heat <- norm.data %>% 
  select(group, test, name, norm, change_to_ref) %>% 
  group_by(test, name, norm) %>% 
  summarise(m=mean(-change_to_ref), .groups = 'drop') %>% 
  mutate(norm=factor(norm, levels=names(norm.colours))) %>%  
  mutate(name=factor(name, levels = rev(c('precision', 'recall', 'AUC')))) %>% 
  mutate(test=factor(test, levels=names(test.colours))) %>% 
  ggplot(aes(x=test, y=name, fill=m)) + 
    geom_tile() + 
    facet_grid(~norm, scale='free', space='free_x') +
    theme_publication() + 
    xlab('') + ylab('') + 
    scale_fill_gradientn(
      colours=rev(c("#193F90", '#3B6FB6', '#8BB8E8', 'white',
                    '#E58F9E', '#D41645', '#A6093D')), 
      limits=c(-0.3, 0.3))

g.norm.box <- norm.data %>% 
  mutate(norm=factor(norm, levels=names(norm.colours))) %>%  
  mutate(name=factor(name, levels = c('precision', 'recall', 'AUC'))) %>% 
  mutate(test=factor(test, levels=names(test.colours))) %>% 
  ggplot(aes(x=test, y=-change_to_ref, fill=test)) + 
    geom_boxplot() +
    facet_grid(name~norm, scale='free_x', space='free_x') +
    scale_fill_manual(values=test.colours) +
    theme_publication(panel.grid = 'major_y') +
    coord_cartesian(ylim=c(-0.5, 0.5))

pdf(here('figures', 'figure_vanilla', 'sim_all_norm.pdf'),
    width = 8, height = 4, useDingbats = FALSE)
print(g.norm.bar)
print(g.norm.box)
print(g.norm.heat)
dev.off()


norm.data %>% 
  group_by(group, subset, test, norm) %>% 
  summarize(overall=sum(change_to_ref), .groups = 'drop') %>% 
  group_by(test, norm, .groups='drop') %>%
  summarise(m=mean(overall)) %>% 
  arrange(m) %>% 
  group_by(test) %>% 
  slice(1)

norm.rank <- .f_find_norms(df.sim.all)

g <- df.sim.all %>% 
  filter(group=='ab4_prev2') %>% 
  mutate(g=paste0(test, '-', norm)) %>% 
  filter(g %in% norm.rank$g) %>% 
  select(subset, test, precision, recall, AUC, fdr) %>% 
  mutate(test=factor(test, levels = names(test.colours))) %>% 
  pivot_longer(c(precision, recall, AUC, fdr)) %>% 
  mutate(value=case_when(name!='fdr'~value, TRUE~log10(value+1e-03))) %>% 
  mutate(lt=test=='metagenomeSeq2') %>% 
  ggplot(aes(x=as.factor(subset), y=value)) + 
    geom_line(aes(group=test, colour=test, lty=lt)) + 
    scale_colour_manual(values=unlist(test.colours)) +
    facet_wrap(~name, scales = 'free_y') + 
    theme_presentation(panel.grid = 'major_y')

# FDR control over all conditions


df.sim.all %>% 
  mutate(g=paste0(test, '-', norm)) %>% 
  filter(g %in% norm.rank$g) %>% 
  filter(!ab %in% c('ab6', 'ab7')) %>%
  filter(subset!=800) %>%
  mutate(x=fdr > 0.05) %>% 
  # ggplot(aes(x=precision, y=recall, col=group)) + 
    # geom_point() + geom_line(aes(group=as.factor(subset))) + 
    # facet_wrap(~test)
  # filter(subset %in% c(50, 100, 200)) %>%
  # group_by(test) %>%
  # summarise(m=sum(x))
  ggplot(aes(x=group, y=as.factor(subset), fill=x)) + 
    geom_tile() + 
    facet_wrap(~test)

  select(subset, test, precision, recall, AUC, fdr) %>% 
  mutate(test=factor(test, levels = names(test.colours))) %>% 
  pivot_longer(c(precision, recall, AUC, fdr)) %>% 
  mutate(value=case_when(name!='fdr'~value, TRUE~log10(value+1e-03))) %>%

        
ranks.test <- df.sim.all %>% 
  filter(!ab %in% c('ab6', 'ab7')) %>% 
  mutate(g=paste0(test, '-', norm)) %>% 
  filter(g %in% norm.rank$g) %>% 
  filter(subset > 30 & subset < 300) %>%
  group_by(type) %>% 
  summarise(m.fdr=mean(fdr), m.auc=mean(AUC), n.fdr=sum(fdr > 0.05), 
            .groups = 'drop') %>% 
  # arrange(n.fdr, desc(m.auc)) %>% 
  arrange(m.fdr) %>% 
  separate(type, into=c('test', 'norm'), sep='-')

g.line <- df.sim.all %>% 
  filter(group=='ab4_prev1') %>% 
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

.f_find_test_ranks <- function(df.res, n){
  df.res %>% 
    mutate(g=paste0(test, '-', norm)) %>% 
    filter(g %in% n$g) %>% 
    filter(subset %in% c(50, 100, 200)) %>%
    filter(!ab %in% c('ab6', 'ab7')) %>% 
    group_by(test) %>% 
    summarise(n.fdr=mean(fdr > 0.1), mean.auroc=mean(AUC)) %>% 
    mutate(n.fdr=as.numeric(n.fdr > 0.1)) %>% 
    arrange(n.fdr, desc(mean.auroc))
}

# other 
df.sim.MMH <- .f_load_and_preprocess('sim_MMH')
df.sim.Weiss <- .f_load_and_preprocess('sim_Weiss')
df.sim.dir <- .f_load_and_preprocess('sim_dirichlet')
df.sim.beta <- .f_load_and_preprocess('sim_betabin')
df.sim.neg <- .f_load_and_preprocess('sim_negbin')
df.sim.spd <- .f_load_and_preprocess('sim_sparseDOSSA')


norm.rank.MMH <- .f_find_norms(df.sim.MMH)
ranks.MMH <- .f_find_test_ranks(df.sim.MMH, norm.rank.MMH) %>% 
  mutate(type='MMH')

norm.rank.Weiss <- .f_find_norms(df.sim.Weiss)
ranks.Weiss <- .f_find_test_ranks(df.sim.Weiss, norm.rank.Weiss) %>% 
  mutate(type='Weiss')

norm.rank.dir <- .f_find_norms(df.sim.dir)
ranks.dir <- .f_find_test_ranks(df.sim.dir, norm.rank.dir) %>% 
  mutate(type='dir')

norm.rank.beta <- .f_find_norms(df.sim.beta)
ranks.beta <- .f_find_test_ranks(df.sim.beta, norm.rank.beta) %>% 
  mutate(type='betabin')

norm.rank.neg <- .f_find_norms(df.sim.neg)
ranks.neg <- .f_find_test_ranks(df.sim.neg, norm.rank.neg) %>% 
  mutate(type='negbin')

norm.rank.spd <- .f_find_norms(df.sim.spd)
ranks.spd <- .f_find_test_ranks(df.sim.spd, norm.rank.spd) %>% 
  mutate(type='sparseDOSSA')

ranks.all <- .f_find_test_ranks(df.sim.all, norm.rank) %>% 
  mutate(type='sim_all')


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
