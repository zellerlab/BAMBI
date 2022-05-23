# ##############################################################################
#
## Simulations with artificial confounding
#
# ##############################################################################

library("here")
library("tidyverse")
library("progress")
library("ggembl")
library("ggpattern")

devtools::load_all(simba.loc)
test.colours <- yaml::yaml.load_file(here('files', 'test_colours.yml'))

# phi coefficient (from the psych packages)
phi <- function(t){  
  # expects: t is a 2 x 2 matrix or a vector of length(4)
  stopifnot(prod(dim(t)) == 4 || length(t) == 4)
  if(is.vector(t)) t <- matrix(t, 2)
  r.sum <- rowSums(t)
  c.sum <- colSums(t)
  total <- sum(r.sum)
  r.sum <- r.sum/total
  c.sum <- c.sum/total
  v <- prod(r.sum, c.sum)
  phi <- (t[1,1]/total - c.sum[1]*r.sum[1]) /sqrt(v)
  names(phi) <- NULL
  return(phi)
}

# ##############################################################################
# Artifical confounder viz

df.plot.all <- list()

sim.files <- list.files(here('simulations', 'conf_sim'), 
                        pattern='artificial.*.h5')
bias.values <- h5read(here('simulations', 'conf_sim', sim.files[1]), 
                      name='/simulation_parameters')
bias.values <- bias.values$conf_bias_test_idx$bias

conf.bar.plot <- list()
temp <- h5read(here('simulations', 'conf_sim', sim.files[3]), 
               name='ab1_prev1_rep1')
for (b in seq_along(bias.values)){
  idx <- temp$conf_bias_test_idx[[b]]$subset_400
  all <- map(seq_len(nrow(idx))[-1], .f = function(x){
    tibble(label=idx[1,], conf=temp$conf_label[idx[x,]]) %>% 
      group_by(label, conf) %>% 
      tally() %>% 
      group_by(label) %>% 
      mutate(n.all=sum(n)) %>% 
      mutate(prop=n/n.all) %>% 
      select(label, conf, prop)
  }) %>% 
    bind_rows() %>% group_by(label, conf) %>% 
    summarise(m=mean(prop), .groups = 'drop') %>% 
    mutate(bias=b, bias.value=bias.values[b])
  conf.bar.plot[[b]] <- all
}
conf.bar.plot <- bind_rows(conf.bar.plot)
conf.bar <- conf.bar.plot %>% 
  ggplot(aes(x=as.factor(label), y=m, fill=as.factor(conf))) + 
    geom_bar(stat='identity') + 
    facet_grid(~bias) +
    theme_publication() + 
    xlab('') + ylab('Proportion')
ggsave(conf.bar, filename = here('figures', 'figure_conf_artificial', 
                                 'conf_bar_plot.pdf'),
       width = 5, height = 3, useDingbats=FALSE)  


fn.df.plot <- here('figures', 'figure_conf_artificial', 'sim_fold_change.tsv')
if (!file.exists(fn.df.plot)){
  for (x in sim.files){
    message(x)
    # number of unique samples
    pb <- progress_bar$new(total=50)
    for (rep in seq_len(50)){
      temp <- h5read(here('simulations', 'conf_sim', x),  
                     name=paste0('ab3_prev2_rep', rep))
      feat <- prop.table(temp$features, 2)
      rownames(feat) <- temp$feature_names
      colnames(feat) <- temp$sample_names
      
      conf.label <- temp$conf_label
      names(conf.label) <- temp$sample_names
      
      for (b in names(temp$conf_bias_test_idx)){
        idx.mat <- temp$conf_bias_test_idx[[b]]$subset_400
        idx.ctr <- idx.mat[2,idx.mat[1,]==1]
        idx.case <- idx.mat[2,idx.mat[1,]==-1]
        
        # label effect
        gfc.label <- vapply(seq_len(nrow(feat)), FUN = function(x){
          q.ctr <- quantile(log10(feat[x,idx.ctr] + 1e-05), 
                            probs=seq(from=0.05, to=0.95, by=0.05))
          q.case <- quantile(log10(feat[x,idx.case] + 1e-05), 
                             probs=seq(from=0.05, to=0.95, by=0.05))
          mean(q.ctr-q.case)
        }, FUN.VALUE = double(1))
        
        # confounder effect
        idx.c.ctr <- idx.mat[2,][which(conf.label[idx.mat[2,]]==1)]
        idx.c.case <- idx.mat[2,][which(conf.label[idx.mat[2,]]==-1)]
        
        gfc.conf <- vapply(seq_len(nrow(feat)), FUN = function(x){
          q.ctr <- quantile(log10(feat[x,idx.c.ctr] + 1e-05), 
                            probs=seq(from=0.05, to=0.95, by=0.05))
          q.case <- quantile(log10(feat[x,idx.c.case] + 1e-05), 
                             probs=seq(from=0.05, to=0.95, by=0.05))
          mean(q.ctr-q.case)
        }, FUN.VALUE = double(1))
        
        # calculate other things
        # number of unique samples
        n.samples <- map(seq_len(nrow(idx.mat))[-1], 
                        .f = function(x){length(unique(idx.mat[x,]))}) %>% 
          unlist() %>% mean
        # phi coefficient between label and confounder
        phi.coef <- map(seq_len(nrow(idx.mat))[-1], 
                        .f=function(x){
                          phi(as.matrix(table(idx.mat[1,],
                                              conf.label[idx.mat[x,]])))
                        }) %>% 
          unlist() %>% mean()
        df.plot.all[[(length(df.plot.all) + 1)]] <- tibble(
          id=rownames(feat), label=gfc.label, conf=gfc.conf) %>% 
          mutate(type=case_when(id %in% temp$marker_idx~'label',
                                id %in% temp$marker_idx_conf~'conf',
                                TRUE~'background')) %>% 
          mutate(n.unique=n.samples, phi=phi.coef) %>% 
          mutate(bias=b) %>% 
          mutate(sim=x) %>% 
          mutate(rep=rep)
      }
      pb$tick()
    }
  }
  df.plot.all <- bind_rows(df.plot.all)
  write_tsv(df.plot.all, file = fn.df.plot)
} else {
  df.plot.all <- read_tsv(fn.df.plot)
}

g.fc <- df.plot.all %>%
  filter(rep==10, bias %in% c('bias_1', 'bias_3', 'bias_6')) %>% 
  mutate(type=factor(type, levels=c('background', 'conf', 'label'))) %>% 
  arrange(type) %>% 
  ggplot(aes(x=label, y=-conf, col=type)) + 
    geom_point() + 
    facet_grid(sim~bias) + 
    xlab("gFC for the label") +
    ylab("gFC for the confounder") +
    theme_publication(panel.grid = 'major') + 
    scale_colour_manual(values=c(alpha('#707372', alpha=0.5), 
                                 '#FFA300', '#A8C700'), 
                        guide=FALSE)

g.cor <- df.plot.all %>% 
  group_by(bias, sim, rep) %>% 
  summarise(correlation=cor(label, -conf), 
            n.unique=unique(n.unique),
            phi=unique(phi), .groups = 'drop') %>% 
  pivot_longer(cols=c(correlation, n.unique, phi)) %>% 
  mutate(value=case_when(name=='n.unique'~value/400, 
                         name=='phi'~-value,TRUE~value)) %>% 
  group_by(bias, sim, name) %>% 
  summarise(m=mean(value), s=sd(value), .groups = 'drop') %>% 
  mutate(bias=str_remove(bias, 'bias_')) %>% 
  mutate(bias.v=bias.values[as.numeric(bias)]) %>%
  ggplot(aes(x=bias.v, y=m)) + 
    geom_line(aes(group=name, col=name)) +
    geom_ribbon(aes(ymin=m-s, ymax=m+s, group=name), alpha=0.2) +
    theme_publication(panel.grid = 'major_y') + 
    xlab('') + ylab('Correlation of gFC values') + 
    geom_abline(slope=2, intercept = -1, colour='darkgrey') + 
    guides(colour='none') + 
    facet_grid(sim~.)

ggsave(g.fc, filename = here('figures', 'figure_conf_artificial', 'gFC.pdf'),
       width = 6, height = 6, useDingbats=FALSE)
ggsave(g.cor, filename = here('figures','figure_conf_artificial', 'cor.pdf'),
       width = 6, height = 7, useDingbats=FALSE)

# ##############################################################################
# load results
fn.final <- here('figures', 'figure_conf_artificial', 
                 paste0('test_results.tsv'))
  
if (file.exists(fn.final)){
  df.res <- read_tsv(fn.final)
} else {
    fn.res <- list.files(here('test_results'), 
                         pattern = paste0('artificial', '.+tsv$'), 
                         recursive = TRUE)
    df.res <- map(fn.res, .f = function(x){
      res <- read_tsv(here('test_results', x)) %>% 
        separate(group, into=c('group', 'repetition'), sep='_rep') %>% 
        separate(subset, into=c('bias', 'subset'), sep='-') %>% 
        mutate(subset=str_remove(subset, 'subset_')) %>% 
        mutate(subset=as.factor(as.numeric(subset))) %>% 
        mutate(PR = 1-FDR) %>% 
        select(-job.id, -time.running) %>% 
        group_by(group, subset, test, norm, bias, problem) %>% 
        summarise(precision=mean(PR), sd.precision=sd(PR),
                  recall=mean(R), sd.recall=sd(R),
                  AUC=mean(auroc), sd.AUC=sd(auroc),
                  fdr=mean(FDR), sd.fdr=sd(FDR),
                  .groups='drop') %>%
        separate(group, into=c('ab','prev'), 
                 sep = '_', remove=FALSE)
        write_tsv(res, file=fn.final, append = file.exists(fn.final))
      })
      df.res <- read_tsv(fn.final)
}


# plot FDR/AUC for confounder-corrected and naive tests according to the
# bias term:
# as line plots
g.lines <- df.res %>% 
  filter(group=='ab3_prev4') %>% 
  filter(subset==200) %>% 
  mutate(corrected=str_detect(test, '_conf')) %>% 
  mutate(test=str_remove(test, '_conf')) %>% 
  select(test, bias, precision, sd.precision, recall, sd.recall, 
         AUC, sd.AUC, corrected, problem) %>% 
  pivot_longer(cols=c(precision, AUC, recall), 
               names_to = 'type', values_to = 'value') %>% 
  pivot_longer(cols=c(sd.precision, sd.AUC, sd.recall), names_to = 'sd.type', 
               values_to = 'sd') %>% 
  mutate(sd.type=str_remove(sd.type, 'sd.')) %>% 
  mutate(test=case_when(test=='lme'~'lm', TRUE~test)) %>% 
  filter(type==sd.type) %>% 
  mutate(bias=str_remove(bias, 'bias_')) %>% 
  mutate(bias.v=bias.values[as.numeric(bias)]) %>% 
  mutate(g=paste0(type, test, corrected)) %>% 
  mutate(value=case_when(type=='AUC'~(value-0.5)/0.5, TRUE~value)) %>%
  ggplot(aes(x=bias.v, y=value)) + 
  geom_line(aes(group=g, col=type, linetype=corrected)) + 
  geom_ribbon(aes(ymax=value+sd, ymin=value-sd, 
                  group=g, fill=type), alpha=0.1) + 
  facet_grid(test~problem) + 
  scale_colour_manual(values=c('#3B6FB6','#D41645', '#18974C')) +
  scale_fill_manual(values=c('#3B6FB6','#D41645', '#18974C')) +
  scale_linetype_manual(values=c(3, 1)) +
  theme_publication()
ggsave(g.lines, filename = here('figures', 'figure_conf_artificial', 
                                'lines.pdf'),
       width = 6, height = 6, useDingbats=FALSE)

# as bar plots for selected values
g.bar <- df.res %>% 
  filter(group=='ab3_prev4') %>% 
  filter(problem=='sim_artificial_conf_0.1') %>% 
  filter(subset==200) %>% 
  mutate(corrected=str_detect(test, '_conf')) %>% 
  mutate(test=str_remove(test, '_conf')) %>% 
  select(test, bias, precision, sd.precision, recall, sd.recall, 
         AUC, sd.AUC, corrected, problem) %>% 
  pivot_longer(cols=c(precision, AUC, recall), 
               names_to = 'type', values_to = 'value') %>% 
  pivot_longer(cols=c(sd.precision, sd.AUC, sd.recall), names_to = 'sd.type', 
               values_to = 'sd') %>% 
  mutate(sd.type=str_remove(sd.type, 'sd.')) %>% 
  mutate(test=case_when(test=='lme'~'lm', TRUE~test)) %>% 
  filter(type==sd.type) %>% 
  mutate(bias=str_remove(bias, 'bias_')) %>% 
  mutate(bias.v=bias.values[as.numeric(bias)]) %>% 
  mutate(g=paste0(type, test, corrected)) %>% 
  mutate(type=factor(type, levels = c('precision', 'recall', 'AUC'))) %>% 
  filter(bias.v %in% c(0.5, 0.7, 0.9)) %>% 
  ggplot(aes(x=test, fill=test, y=value, alpha=corrected)) +
    geom_bar(stat='identity', position = position_dodge()) +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd),
                  position = position_dodge(width = 0.9), width=0.2) +
    facet_grid(type~as.factor(bias.v)) +
    theme_publication(panel.grid = 'major_y') +
    scale_alpha_manual(values=c(1, 1)) +
    scale_fill_manual(values=test.colours) +
    geom_hline(yintercept = 0.9, colour='black')
ggsave(g.bar, filename = here('figures', 'figure_conf_artificial', 
                              'barplot.pdf'),
       width = 6, height = 6, useDingbats=FALSE)


# heatmaps maybe?
group.colour <- c(
  "#fef5a7", "#fdf07c", "#fdeb50", "#FDE725FF",
  "#d2efb4", "#bbe78e", "#a5df69", "#8FD744FF",
  "#aee2c9", "#85d3ae", "#5dc593", "#35B779FF",
  "#a6d2d1", "#79bcba", "#4da6a3", "#21908CFF",
  "#acc2d1", "#83a4bb", "#5a86a4", "#31688EFF",
  "#b4b0cd", "#8e88b4", "#69619b", "#443a83FF",
  "#b499ba", "#8e6698", "#693376", "#440154FF"
)
df.res %>% 
  mutate(corrected=str_detect(test, '_conf')) %>% 
  mutate(test=str_remove(test, '_conf')) %>%  
  filter(bias %in% c('bias_1','bias_3','bias_6')) %>% 
  filter(problem == 'sim_artificial_conf_0.2') %>% 
  filter(subset=='200') %>%
  mutate(test=case_when(test=='lme'~'lm', TRUE~test)) %>% 
  arrange(corrected) %>% 
  ggplot(aes(x=precision, y=recall)) + 
    geom_point(aes(fill=corrected), pch=21) +
    geom_line(aes(group=group, col=group))+
    facet_grid(test~bias) +
    scale_x_continuous(limits=c(0,1)) + 
    scale_y_continuous(limits=c(0,1)) + 
    scale_colour_manual(values=group.colour) + 
    theme_publication(panel.grid = 'major') + 
    scale_fill_manual(values=c('red', 'black'))
