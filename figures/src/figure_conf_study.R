# ##############################################################################
#
## Simulations with study confounding
#
# ##############################################################################

library("here")
library("tidyverse")
library("progress")
library("ggembl")
library("cowplot")

devtools::load_all(simba.loc)
test.colours <- yaml::yaml.load_file(here('files', 'test_colours.yml'))
dataset.shapes <- c('Zeevi_WGS'=1, 'Schirmer_WGS'=3, 'TwinsUK_WGS'=5)

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

# find all study confounder simulations
sim.files <- list.files(here('simulations', 'conf_sim'), 
                        pattern='sim_study_conf')
for (i in sim.files){
  tag <- str_remove(i, 'sim_study_conf_') %>% 
    str_remove('.h5')
  message(tag)
  if (!dir.exists(here('figures', 'figure_conf_study', tag))){
    dir.create(here('figures', 'figure_conf_study', tag))
  }
  sim.params <- h5read(here('simulations', 'conf_sim', i),
                       name='/simulation_parameters')
  bias.values <- sim.params$conf_bias_test_idx$bias
  
  # ############################################################################
  # gFC and PCoA plots
  df.plot.gfc <- list()
  
  temp <- h5read(here('simulations', 'conf_sim', i),
                 name='ab3_prev3_rep10')
  feat <- temp$features
  colnames(feat) <- temp$sample_names
  rownames(feat) <- temp$feature_names
  feat.rel <- prop.table(feat, 2)
  pco.res <- labdsv::pco(vegan::vegdist(t(feat)))
  colnames(pco.res$points) <- c('PCo1', 'PCo2')
  axis.text <- sprintf(fmt='%.2f', 
                       pco.res$eig[1:2]/sum(pco.res$eig[
                         pco.res$eig > 0]) * 100)
  df.pcoa <- pco.res$points %>% 
    as_tibble(rownames = 'Sample_ID') %>% 
    mutate(conf=temp$conf_label)
  g.pcoa <- list()
  g.phis <- tibble(bias=double(0), phi=double(0))
  message("Calculating PCoAs and gFC values across the biases")
  pb <- progress_bar$new(total=length(bias.values))
  for (i.b in seq_along(bias.values)){
    # index
    idx <- temp$conf_bias_test_idx[[i.b]]$subset_400
    
    # calculate phis
    g.phis <- bind_rows(g.phis, 
                        tibble(bias=i.b, phi=unlist(map(
                          seq_len(nrow(idx))[-1], 
                          .f=function(x){
                            -phi(table(idx[1,], temp$conf_label[idx[x,]]))}))))
  
    # make a PCoA
    g <- df.pcoa %>% 
      mutate(bias=i.b) %>% 
      as.data.frame()
    rownames(g) <- g$Sample_ID
    g <- g[temp$sample_names[idx[13,]],]
    g$label <- idx[1,]
    g.pcoa[[i.b]] <- g
    
    # calculate gFCs
    gfc <- map(seq_len(nrow(idx))[-1], .f=function(rep){
      vapply(rownames(feat.rel), FUN=function(x){
        f <- log10(feat.rel[x,idx[rep,]] + 1e-05)
        l <- idx[1,]
        c <- temp$conf_label[idx[rep,]]
        q.p <- quantile(f[l==1], probs = seq(.1, .9, .05))
        q.n <- quantile(f[l==-1], probs = seq(.1, .9, .05))
        q.cp <- quantile(f[c==1], probs = seq(.1, .9, .05))
        q.cn <- quantile(f[c==2], probs = seq(.1, .9, .05))
        l <- mean(q.p-q.n)
        c <- mean(q.cp-q.cn)
        return(c('label'=l, 'conf'=c))}, FUN.VALUE = double(2)) %>%
        t() %>%
        as_tibble() %>%
        mutate(motu=rownames(feat.rel), bias=i.b) %>%
        mutate(type=motu %in% temp$marker_idx) %>%
        mutate(repetition=rep)}) %>%
      bind_rows()
    df.plot.gfc[[i.b]] <- gfc
    pb$tick()
  }

  # phis
  g.phis %>% group_by(bias) %>% summarise(m=mean(phi), s=sd(phi)) %>% print

  # PCoA
  g.pcoa <- g.pcoa %>% 
    bind_rows()
  datasets <- str_split(tag, '-')[[1]]
  g1 <- g.pcoa %>% 
    mutate(conf.dataset=datasets[conf]) %>% 
    ggplot(aes(x=PCo1, y=PCo2, col=as.factor(label))) + 
    geom_point(aes(shape=conf.dataset)) + 
    theme_publication() + 
    facet_grid(~bias) +
    scale_colour_manual(values=c('#734595', '#F49E17')) + 
    scale_shape_manual(values=dataset.shapes) +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    guides(colour='none') +
    xlab(paste0('PCo1 [', axis.text[1], '%]')) + 
    ylab(paste0('PCo2 [', axis.text[2], '%]'))
  g2 <- g.pcoa %>% 
    mutate(conf.dataset=datasets[conf]) %>% 
    ggplot(aes(x=as.factor(label), y=PCo1)) + 
    geom_boxplot(aes(fill=as.factor(label))) +
    geom_jitter(aes(col=as.factor(label), shape=conf.dataset), width = 0.1) +
    theme_publication() +
    facet_grid(~bias) + 
    scale_fill_manual(values = alpha(c('#734595', '#F49E17'), alpha=0.5)) + 
    scale_colour_manual(values=c('#734595', '#F49E17')) + 
    scale_shape_manual(values=dataset.shapes) +
    guides(colour='none', shape='none', fill='none') + 
    coord_flip() + 
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    xlab('') + ylab('')
  g <- plot_grid(g1, g2, rel_heights = c(0.8, 0.2), ncol = 1)
  ggsave(g, filename = here('figures', 'figure_conf_study', 
                            tag, 'pcoa.pdf'),
       width = 10, height = 4, useDingbats=FALSE)

  # gFC
  df.plot.gfc <- bind_rows(df.plot.gfc)

  # correlation between gFC values
  df.plot.gfc %>% 
    group_by(bias, repetition) %>% 
    summarise(m=cor(abs(label), abs(conf)), .groups = 'drop_last') %>% 
    summarise(x=mean(m), s=sd(m))

  g <- df.plot.gfc %>% 
    filter(repetition==30) %>% 
    arrange(type) %>% 
    ggplot(aes(x=abs(label), y=abs(conf), col=type)) + 
    geom_point() +
    facet_grid(~bias) + 
    theme_publication() + 
    geom_abline(slope = 1, intercept = 0) +
    scale_colour_manual(values=c(alpha('#707372', alpha=0.5), '#A8C700'), 
                        guide=FALSE)
  ggsave(g, filename = here('figures', 'figure_conf_study', tag, 'gFC.pdf'),
         width = 10, height = 4, useDingbats=FALSE)

  # ############################################################################
  # test results
  fn.test <- here('figures', 'figure_conf_study', 
                  tag, 'test_results.tsv')
  if (!file.exists(fn.test)){
    
    test.files <- list.files(here('test_results'),
                             pattern = paste0('sim_study_conf_', tag), 
                             recursive = TRUE)
    if (length(test.files) == 0){
      message("Tests have not yet been completed for simuation: ", tag)
      next()
    }
    message("+ reading and combining test results")
    df.test <- map(test.files,.f=function(x){
      read_tsv(here('test_results', x),  col_types = cols()) %>% 
        select(-job.id, -time.running, -problem, -mem.used) %>% 
        separate(group, into=c('group', 'repetition'), sep='_rep') %>% 
        separate(subset, into=c('bias', 'subset'), sep='-') %>% 
        mutate(PR = 1-FDR) %>% 
        group_by(group, subset, bias, test, norm, adjust) %>% 
        summarise(precision=mean(PR), sd.precision=sd(PR),
                  recall=mean(R), sd.recall=sd(R),
                  AUC=mean(auroc), sd.AUC=sd(auroc),
                  fdr=mean(FDR), sd.fdr=sd(FDR),
                  .groups='drop') %>% 
        separate(group, into=c('ab','prev'), sep = '_', remove=FALSE)}) %>% 
      bind_rows()
    message("+ finished reading and combining test results")
    write_tsv(df.test, file = fn.test)
    df.comb <- df.test
    rm(df.test)
  } else {
    df.comb <- read_tsv(file = fn.test, col_types = cols())
  }

  
  # single effect size plot?
  g.lines.nice <- df.comb %>% 
    filter(subset=='subset_200', group=='ab4_prev2', adjust=='BH') %>% 
    mutate(corrected=str_detect(test, '_conf')) %>% 
    mutate(test=str_remove(test, '_conf')) %>% 
    select(group, test, bias, precision, sd.precision, recall, sd.recall, 
           AUC, sd.AUC, corrected) %>% 
    pivot_longer(cols=c(precision, AUC, recall), 
                 names_to = 'type', values_to = 'value') %>% 
    pivot_longer(cols=c(sd.precision, sd.AUC, sd.recall), names_to = 'sd.type', 
                 values_to = 'sd') %>% 
    mutate(sd.type=str_remove(sd.type, 'sd.')) %>% 
    mutate(test=case_when(test=='lme'~'lm', TRUE~test)) %>% 
    filter(type==sd.type) %>% 
    mutate(bias=as.numeric(str_remove(bias, 'bias_'))) %>% 
    mutate(g=paste0(test, corrected, type)) %>% 
    mutate(value=case_when(type=='AUC'~(value-0.5)/0.5, TRUE~value)) %>% 
    mutate(bias.v=bias.values[bias]) %>% 
    mutate(type=factor(type, levels = c('precision', 'recall', 'AUC'))) %>% 
    mutate(bias.v=case_when(test=='limma'~bias.v-0.01, 
                            test=='wilcoxon'~bias.v+0.01,
                            TRUE~bias.v)) %>% 
    ggplot(aes(x=bias.v, y=value)) + 
    geom_point(aes(col=test)) + 
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd, col=test, width=0)) +
    geom_line(aes(group=g, col=test)) +
    facet_grid(corrected~type) + 
    scale_colour_manual(values=unlist(test.colours)) + 
    theme_publication(panel.grid = 'major') +
    scale_y_continuous(limits=c(0,1)) +
    ylab('') + xlab('Confounding strength') +
    NULL
  ggsave(g.lines.nice, filename = here('figures', 'figure_conf_study', 
                                       tag, 'lines_nice.pdf'),
         width = 7, height = 4, useDingbats=FALSE)    
    
  
  
  # plot line plots (for the supplement)
  # being dependent on effect size plots for supplement?
  g.lines <- df.comb %>% 
    filter(subset=="subset_200") %>% 
    mutate(corrected=str_detect(test, '_conf')) %>% 
    mutate(test=str_remove(test, '_conf')) %>% 
    filter(prev=='prev2') %>%
    select(group, test, bias, precision, sd.precision, recall, sd.recall, 
           AUC, sd.AUC, corrected, adjust) %>% 
    pivot_longer(cols=c(precision, AUC, recall), 
                 names_to = 'type', values_to = 'value') %>% 
    pivot_longer(cols=c(sd.precision, sd.AUC, sd.recall), names_to = 'sd.type', 
                 values_to = 'sd') %>% 
    mutate(sd.type=str_remove(sd.type, 'sd.')) %>% 
    mutate(test=case_when(test=='lme'~'lm', TRUE~test)) %>% 
    filter(type==sd.type) %>% 
    mutate(bias=as.numeric(str_remove(bias, 'bias_'))) %>% 
    mutate(g=paste0(test, corrected, type)) %>% 
    mutate(value=case_when(type=='AUC'~(value-0.5)/0.5, TRUE~value)) %>% 
    mutate(bias.v=bias.values[bias]) %>% 
    # mutate(sd=case_when(type=='AUC'~(sd-0.5)/0.5, TRUE~sd)) %>%
    ggplot(aes(x=bias.v, y=value)) + 
      geom_line(aes(group=g, col=type, linetype=corrected)) + 
      geom_ribbon(aes(ymin=value-sd, ymax=value+sd, group=g, fill=type), 
                  alpha=0.3) +
      facet_grid(group~adjust+test) + 
      theme_publication() +
      coord_cartesian(ylim=c(0,1)) +
      scale_colour_manual(values=c('#3B6FB6','#D41645', '#18974C')) +
      scale_fill_manual(values=c('#3B6FB6','#D41645', '#18974C')) +
      scale_linetype_manual(values=c(3, 1)) +
      NULL
  # ggsave(g.lines, filename = here('figures', 'figure_conf_study', 
  #                                'lines.pdf'),
  #       width = 10, height = 8, useDingbats=FALSE)    

  # plot bar plots (for the main text)
  g <- df.comb %>% 
    filter(subset=="subset_200") %>% 
    mutate(corrected=str_detect(test, '_conf')) %>% 
    mutate(test=str_remove(test, '_conf')) %>% 
    filter(group=='ab4_prev2') %>% 
    select(test, bias, precision, sd.precision, recall, sd.recall, 
           AUC, sd.AUC, corrected, adjust) %>% 
    pivot_longer(cols=c(precision, AUC, recall), 
                 names_to = 'type', values_to = 'value') %>% 
    pivot_longer(cols=c(sd.precision, sd.AUC, sd.recall), names_to = 'sd.type', 
                 values_to = 'sd') %>% 
    mutate(sd.type=str_remove(sd.type, 'sd.')) %>% 
    mutate(test=case_when(test=='lme'~'lm', TRUE~test)) %>% 
    filter(type==sd.type) %>% 
    mutate(bias=str_remove(bias, 'bias_')) %>% 
    ggplot(aes(x=test, y=value, alpha=corrected, fill=test)) + 
      geom_bar(stat='identity', position = position_dodge()) +
      geom_errorbar(aes(ymin=value-sd, ymax=value+sd),
                    position = position_dodge(width = 0.9), width=0.2) +
      facet_grid(type~adjust+as.factor(bias)) +
      theme_publication(panel.grid = 'major_y') +
      scale_alpha_manual(values=c(1, 1)) +
      scale_fill_manual(values=test.colours) +
      geom_hline(yintercept = 0.9, colour='black')
  # ggsave(g, filename = here('figures', 'figure_conf_study', 
  #                          'barplots.pdf'),
  #       width = 10, height = 5, useDingbats=FALSE)
}
