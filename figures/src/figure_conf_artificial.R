# ##############################################################################
#
## Simulations with artificial confounding
#
# ##############################################################################

library("here")
library("tidyverse")
library("progress")
library("ggembl")

set.seed(2022)

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

sim.files <- list.files(here('simulations', 'conf_sim'), 
                        pattern='artificial.*.h5')
df.test.all <- list()
for (fn.sim in sim.files){
  tag <-  str_remove(fn.sim, '.h5')
  message(tag)
  if (!dir.exists(here('figures', 'figure_conf_artificial', tag))){
    dir.create(here('figures', 'figure_conf_artificial', tag))
  }
  
  bias.values <- h5read(here('simulations', 'conf_sim', fn.sim), 
                        name='/simulation_parameters')
  bias.values <- bias.values$conf_bias_test_idx$bias

  conf.bar.plot <- list()
  temp <- h5read(here('simulations', 'conf_sim', fn.sim), 
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
      xlab('') + ylab('Proportion') +
     scale_fill_manual(values=c('#00589C', '#5EADCE'))
  ggsave(conf.bar, filename = here('figures', 'figure_conf_artificial', 
                                   tag, 'conf_bar_plot.pdf'),
         width = 5, height = 3, useDingbats=FALSE)  


  fn.df.plot <- here('figures', 'figure_conf_artificial', 
                     tag, 'sim_fold_change.tsv')
  if (!file.exists(fn.df.plot)){
    # number of unique samples
    df.plot.all <- list()
    pb <- progress_bar$new(total=50)
    for (rep in seq_len(50)){
      temp <- h5read(here('simulations', 'conf_sim', fn.sim),  
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
          mutate(sim=str_remove(fn.sim, '.h5')) %>% 
          mutate(rep=rep)
      }
      pb$tick()
    }
    df.plot.all <- bind_rows(df.plot.all)
    write_tsv(df.plot.all, file = fn.df.plot)
  } else {
    df.plot.all <- read_tsv(fn.df.plot, col_types = cols())
  }

  # plot some examples!
  # if (FALSE){
    # confounded, no label effect
    temp <- h5read(here('simulations', 'conf_sim', fn.sim),  
                   name='ab3_prev2_rep10')
    feat <- prop.table(temp$features, 2)
    rownames(feat) <- temp$feature_names
    colnames(feat) <- temp$sample_names
    temp.info <- df.plot.all %>% 
      filter(rep==10) %>% 
      left_join(enframe(rowMeans(feat!=0), value='prev', name='id'), by='id')
    
    x1.1 <- temp.info %>% 
      filter(abs(label) < 0.1) %>%
      filter(phi > -0.3) %>% 
      filter(type=='conf') %>% 
      filter(bias=='bias_1') %>% 
      filter(prev > 0.5) %>% 
      arrange(desc(abs(conf))) %>% 
      slice(min(10, nrow(.))) %>% 
      sample_n(1)
    # confounded, also label effect
    x1.2 <- temp.info %>% 
      filter(abs(label) > 0.5) %>% 
      filter(phi < -0.3) %>% 
      filter(type=='conf') %>% 
      filter(bias=='bias_5') %>% 
      arrange(desc(abs(conf))) %>% 
      filter(prev > 0.5) %>%
      slice(min(10, nrow(.))) %>% 
      sample_n(1)
    # not confounded, label effect
    x2.1 <- temp.info %>% 
      filter(abs(conf) < 0.1) %>% 
      filter(phi > -0.3) %>% 
      filter(type=='label') %>% 
      filter(bias=='bias_1') %>% 
      filter(prev > 0.5) %>%
      arrange(desc(abs(label))) %>% 
      slice(min(10, nrow(.))) %>% 
      sample_n(1)
    # label effect, also confounded
    x2.2 <- temp.info %>% 
      filter(abs(conf) > 0.5) %>% 
      filter(phi < -0.3) %>% 
      filter(type=='label') %>% 
      filter(bias=='bias_5') %>% 
      filter(prev > 0.5) %>%
      arrange(desc(abs(label))) %>% 
      slice(min(10, nrow(.))) %>% 
      sample_n(1)
    ex.selected <- bind_rows(x1.1, x1.2, x2.1, x2.2)

    conf.label <- temp$conf_label
    names(conf.label) <- temp$sample_names
    
    df.tmp <- list()
    
    for (a in seq_len(nrow(ex.selected))){
      idx.mat <- temp$conf_bias_test_idx[[ex.selected$bias[a]]]$subset_400
      idx.ctr <- idx.mat[2,idx.mat[1,]==1]
      idx.case <- idx.mat[2,idx.mat[1,]==-1]
      
      df.tmp[[a]] <- tibble(ab=log10(feat[ex.selected$id[a],
                                          c(idx.ctr, idx.case)] + 1e-5),
                            label=c(rep(1, length(idx.ctr)), 
                                    rep(-1, length(idx.case))),
                            conf=conf.label[c(idx.ctr, idx.case)],
                            id=ex.selected$id[a])
    }
    
    
    bias_1 <- idx.mat <- temp$conf_bias_test_idx[['bias_1']]$subset_400
    bias_5 <- idx.mat <- temp$conf_bias_test_idx[['bias_5']]$subset_400
    
    df.tmp <- bind_rows(df.tmp) %>% left_join(ex.selected, by='id')
    
    g.boxes <- df.tmp %>% 
      ggplot(aes(x=as.factor(label.x), y=ab)) + 
        geom_quantile_box(aes(fill=as.factor(label.x))) +
        geom_boxplot(aes(col=as.factor(conf.x)), fill=NA, outlier.shape = NA) +
        geom_jitter(aes(col=as.factor(conf.x)), alpha=0.5, 
                    position = position_jitterdodge()) +
        facet_grid(~id+type+bias, scale='free') + 
        theme_publication(panel.grid='major_y') + 
        scale_fill_manual(values=c('#FFAC1C', '#8159A0')) + 
        scale_colour_manual(values=c('#00589C', '#5EADCE'))
    ggsave(g.boxes, filename = here('figures', 'figure_conf_artificial', tag,
                                           'boxes.pdf'),
           width = 8, height = 4, useDingbats=FALSE)
  if (tag == 'sim_artificial_conf_0.2_Zeevi_WGS') {
    message("+ compute P-values for examples")
    df.p.val <- tibble(id=character(0), test=character(0), 
                       corrected=logical(0), pval=double(0), bias=character(0))
    for (b in c('bias_1', 'bias_5')){
      idx.mat <- temp$conf_bias_test_idx[[b]]$subset_400
      ids <- idx.mat[2,]
      l <- idx.mat[1,]
      names(l) <- colnames(feat)[ids]
      conf <- data.frame(conf=(temp$conf_label + 1)/2 + 1)
      rownames(conf) <- colnames(feat)
      feat.log.TSS <- log10(feat + 1e-05)
      # limma
      df.limma <- test.via.limma(feat.log.TSS[,ids], label = l, conf=NULL)
      df.limma.conf <- test.via.limma(feat.log.TSS[,ids], label=l, conf=conf)
      df.p.val <- df.p.val %>% 
        bind_rows(enframe(df.limma, name='id', value='pval') %>% 
                    mutate(test='limma', corrected=FALSE, bias=b)) %>% 
        bind_rows(enframe(df.limma.conf, name='id', value='pval') %>% 
                    mutate(test='limma', corrected=TRUE, bias=b))
      # lm
      df.lm <- test.lm(feat.log.TSS[ex.selected$id,ids], label = l, conf=NULL)
      df.lm.conf <- test.lme(feat.log.TSS[ex.selected$id,ids], label=l, conf=conf)
      df.p.val <- df.p.val %>% 
        bind_rows(enframe(df.lm, name='id', value='pval') %>% 
                    mutate(test='lm', corrected=FALSE, bias=b)) %>% 
        bind_rows(enframe(df.lm.conf, name='id', value='pval') %>% 
                    mutate(test='lm', corrected=TRUE, bias=b))
      
      # wilcoxon
      df.wilcoxon <- test.wilcoxon(feat[ex.selected$id,ids], label = l, conf=NULL)
      df.wilcoxon.conf <- test.wilcoxon(feat[ex.selected$id,ids], label=l, conf=conf)
      df.p.val <- df.p.val %>% 
        bind_rows(enframe(df.wilcoxon, name='id', value='pval') %>% 
                    mutate(test='wilcoxon', corrected=FALSE, bias=b)) %>% 
        bind_rows(enframe(df.wilcoxon.conf, name='id', value='pval') %>% 
                    mutate(test='wilcoxon', corrected=TRUE, bias=b))
      
      # fastANCOM
      df.fANCOM <- test.fastANCOM(feat[,ids], label = l, conf=NULL)
      df.fANCOM.conf <- test.fastANCOM(feat[,ids], label=l, conf=conf)
      df.p.val <- df.p.val %>% 
        bind_rows(enframe(df.fANCOM, name='id', value='pval') %>% 
                    filter(id %in% ex.selected$id) %>% 
                    mutate(test='fastANCOM', corrected=FALSE, bias=b)) %>% 
        bind_rows(enframe(df.fANCOM.conf, name='id', value='pval') %>% 
                    filter(id %in% ex.selected$id) %>% 
                    mutate(test='fastANCOM', corrected=TRUE, bias=b))
    }
    g.p.val <- df.p.val %>% 
      filter(id %in% ex.selected$id) %>% 
      mutate(label=case_when(pval < 0.05 ~ sprintf(fmt='%.E', pval),
                             TRUE ~ sprintf(fmt='%.2f', pval))) %>% 
      ggplot(aes(x=id, y=test, fill=test)) + 
        geom_tile() + 
        scale_fill_manual(values=unlist(test.colours)) + 
        geom_text(aes(label=label), colour='white') + 
        facet_grid(corrected~bias) + 
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
      xlab('') + ylab('')
    ggsave(g.p.val, filename = here('figures', 'figure_conf_artificial', tag, 
                                    'p_values.pdf'),
           width = 8, height = 6, useDingbats=FALSE)
  }
    
  g.fc <- df.plot.all %>%
    filter(rep==10) %>% 
    mutate(type=factor(type, levels=c('background', 'conf', 'label'))) %>% 
    mutate(id2=paste0(id, '-', bias)) %>% 
    mutate(highlihgt=id2 %in% c(ex.selected %>%
                                  mutate(id2=paste0(id, '-', bias)) %>%
                                  pull(id2))) %>%
    arrange(highlihgt, type) %>%
    ggplot(aes(x=label, y=-conf, fill=type)) + 
      geom_point(aes(col=highlihgt), pch=21) +
      facet_grid(~bias) + 
      xlab("gFC for the label") +
      ylab("gFC for the confounder") +
      theme_publication(panel.grid = 'major') + 
      scale_fill_manual(values=c(alpha('#707372', alpha=0.5), 
                                 '#FFA300', '#A8C700')) + 
      scale_colour_manual(values=c('#FFFFFF00', "black"))

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
      geom_abline(slope=2, intercept = -1, colour='darkgrey') +
      geom_line(aes(group=name, col=name)) +
      geom_ribbon(aes(ymin=m-s, ymax=m+s, group=name), alpha=0.2) +
      theme_publication(panel.grid = 'major_y') + 
      xlab('') + ylab('Correlation of gFC values')

  ggsave(g.fc, filename = here('figures', 'figure_conf_artificial',
                               tag, 'gFC.pdf'),
         width = 9, height = 3, useDingbats=FALSE)
  ggsave(g.cor, filename = here('figures','figure_conf_artificial', 
                                tag, 'cor.pdf'),
         width = 4, height = 2.5, useDingbats=FALSE)
  

  # ############################################################################
  # load results
  fn.final <- here('figures', 'figure_conf_artificial', tag, 'test_results.tsv')
  
  if (file.exists(fn.final)){
    df.res <- read_tsv(fn.final)
  } else {
      fn.res <- list.files(here('test_results'), 
                           pattern = tag,
                           recursive = TRUE)
      if (length(fn.res)==0){
        message("Tests have not yet been finished for simulation: ", tag)
        next()
      }
      df.res <- map(fn.res, .f = function(x){
        res <- read_tsv(here('test_results', x), col_types=cols()) %>% 
          separate(group, into=c('group', 'repetition'), sep='_rep') %>% 
          separate(subset, into=c('bias', 'subset'), sep='-') %>% 
          mutate(subset=str_remove(subset, 'subset_')) %>% 
          mutate(subset=as.factor(as.numeric(subset))) %>% 
          mutate(PR = 1-FDR) %>% 
          select(-job.id, -time.running) %>% 
          group_by(group, subset, test, norm, bias, problem, adjust) %>% 
          summarise(precision=mean(PR), sd.precision=sd(PR),
                    recall=mean(R), sd.recall=sd(R),
                    AUC=mean(auroc), sd.AUC=sd(auroc),
                    fdr=mean(FDR), sd.fdr=sd(FDR),
                    .groups='drop') %>%
          separate(group, into=c('ab','prev'), 
                   sep = '_', remove=FALSE)
          write_tsv(res, file=fn.final, append = file.exists(fn.final))
        })
        df.res <- read_tsv(fn.final, col_types = cols())
  }

  
  # plot FDR/AUC for confounder-corrected and naive tests according to the
  # bias term:
  # as line plots with dots
  
  # important edit!
  # use only those sample sizes that make sense for each dataset (stupid!), 
  # since smaller datasets can be up-sampled with the confounder indices,
  # leading to lots of repeated measurements
  if (str_detect(tag, 'Zeevi_WGS')){
    sample.size=200
  } else {
    sample.size=100
  }
  
  g.lines <- df.res %>% 
    filter(group=='ab4_prev2') %>% 
    filter(subset==sample.size) %>% 
    filter(adjust=='BH') %>%
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
    left_join(df.plot.all %>% group_by(bias) %>% 
                summarize(phi=mean(phi)), by='bias') %>% 
    mutate(bias=str_remove(bias, 'bias_')) %>% 
    mutate(bias.v=bias.values[as.numeric(bias)]) %>% 
    mutate(g=paste0(type, test, corrected)) %>% 
    mutate(value=case_when(type=='AUC'~(value-0.5)/0.5, TRUE~value)) %>%
    mutate(bias.v=case_when(test=='limma'~bias.v-0.01,
                            test=='wilcoxon'~bias.v+0.01,
                            TRUE~bias.v)) %>% 
    mutate(phi=case_when(test=='limma'~phi-0.01,
                         test=='wilcoxon'~phi+0.01,
                         TRUE~phi)) %>% 
    ggplot(aes(x=-phi, y=value, col=test)) + 
    geom_point() + 
    geom_line(aes(group=g, linetype=corrected)) + 
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=0) +
    facet_grid(corrected~type) +
    scale_colour_manual(values=unlist(test.colours)) + 
    theme_publication(panel.grid = 'major_y') + 
    coord_cartesian(ylim=c(0,1)) + 
    scale_linetype_manual(values=c(2,1)) +
    xlab('Strength of confounding') + ylab('')
    
  ggsave(g.lines, filename = here('figures', 'figure_conf_artificial', 
                                  tag, 'lines.pdf'),
         width = 6, height = 4, useDingbats=FALSE)
  
  # over all subset sizes/effect sizes
  # df.res %>% 
  #   filter(adjust=='BH') %>% 
  #   filter(subset==100) %>% 
  #   # mutate(corrected=str_detect(test, '_conf')) %>% 
  #   # mutate(test=str_remove(test, '_conf')) %>% 
  #   # mutate(test=case_when(test=='lme'~'lm', TRUE~test)) %>% 
  #   mutate(precision=precision > 0.9) %>% 
  #   ggplot(aes(x=bias, y=test, fill=precision)) + 
  #     geom_tile() + 
  #     facet_grid(prev~ab)
  # 
  
  
  
  
  # as bar plots for selected values
  g.bar <- df.res %>% 
    filter(group=='ab3_prev4') %>% 
    filter(subset==sample.size) %>% 
    filter(adjust=='BH') %>% 
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
  # ggsave(g.bar, filename = here('figures', 'figure_conf_artificial', 
                                # tag, 'bars.pdf'),
         # width = 6, height = 6, useDingbats=FALSE)
  df.test.all[[tag]] <- df.res %>% 
    mutate(dataset=tag)
}

# ##############################################################################
# plot all results into a single supplementary figure

dataset.colour <- c('0.05_Schirmer_WGS'='#00bddf', 
                    '0.1_Schirmer_WGS'='#0092ac', 
                    '0.2_Schirmer_WGS'='#006679',
                    '0.05_TwinsUK_WGS'='#ffa737', 
                    '0.1_TwinsUK_WGS'='#ff9103', 
                    '0.2_TwinsUK_WGS'='#d07500',
                    '0.05_Zeevi_WGS'='#00e8b0', 
                    '0.1_Zeevi_WGS'='#00b589', 
                    '0.2_Zeevi_WGS'='#008263')

g.all <- df.test.all %>% 
  bind_rows() %>% 
  filter(group=='ab4_prev2') %>% 
  select(subset, test, bias, dataset, adjust, precision, recall, AUC) %>% 
  pivot_longer(c(precision, recall, AUC)) %>% 
  mutate(type=paste0(adjust, '-', name)) %>% 
  filter(type %in% c('BH-precision', 'BH-recall', 'none-AUC')) %>% 
  select(-type) %>% 
  mutate(corrected=str_detect(test, '_conf')) %>% 
  mutate(test=str_remove(test, '_conf')) %>% 
  mutate(test=case_when(test=='lme'~'lm', TRUE~test)) %>% 
  filter(subset==100) %>% 
  mutate(dataset=str_remove(dataset, 'sim_artificial_conf_')) %>% 
  mutate(g=paste0(corrected, dataset)) %>% 
  mutate(bias=str_remove(bias, 'bias_')) %>% 
  mutate(bias.v=bias.values[as.numeric(bias)]) %>% 
  mutate(name=factor(name, levels = c('precision', 'recall', 'AUC'))) %>% 
  ggplot(aes(x=bias.v, y=value, col=dataset)) + 
    geom_line(aes(group=g, linetype=corrected)) +
    geom_point() + 
    facet_grid(name~test) + 
    theme_publication(panel.grid = 'major_y') + 
    xlab('Confounder strength') + ylab('') + 
    ylim(0,1) +
    scale_linetype_manual(values=c(2,1)) + 
    scale_colour_manual(values=dataset.colour)
  
ggsave(g.all, filename = here('figures', 'figure_conf_artificial', 
                              'all_datasets.pdf'),
       width = 170, height = 100, useDingbats=FALSE, unit='mm')



# ##############################################################################
# Example from Sofia's paper?

fn.meta <- '../temp/metadata/meta_Forslund.tsv'
meta <- read.table(fn.meta, sep='\t', stringsAsFactors = FALSE,
                   check.names = FALSE, quote = '')
meta <- meta[meta$GROUP!='IGT',]
meta$METFORMIN[meta$METFORMIN!='Metf'] <- 'NoMetf'

# features
feat <- list()
for (i in c('Karlsson_2013', 'Qin_2012', 'metaHIT')){
  fn.feat <- paste0('../temp/motus/', i, '_motus.tsv')
  x <- read.table(fn.feat, sep='\t', check.names = FALSE,
                  stringsAsFactors = FALSE, quote = '',
                  comment.char = '',
                  header=TRUE) %>% 
    as.matrix()
  # adjust names for the CN cases
  if (i=='Qin_2012') {
    colnames(x) <- str_remove(colnames(x),'bgi-')
  }
  
  feat[[i]] <- x
}

feat <- do.call(cbind, feat)
feat <- feat[,colSums(feat) > 100]

feat <- prop.table(feat, 2)
print(dim(meta))
length(intersect(rownames(meta), colnames(feat)))
meta %>% 
  as_tibble(rownames = 'Sample_ID') %>% 
  filter(Sample_ID %in% colnames(feat)) %>% 
  group_by(COUNTRY) %>% 
  tally()
x <- c('Proteobacteria sp. [ref_mOTU_v25_00095]', 
       'Intestinibacter bartlettii [ref_mOTU_v25_03687]')
g <- as_tibble(t(feat[x,rownames(meta)]), rownames='name') %>% 
  pivot_longer(-name, names_to='species') %>% 
  mutate(value=log10(value+1e-05)) %>% 
  full_join(as_tibble(meta, rownames='name'), by='name') %>% 
  ggplot(aes(x=GROUP, y=value)) +
    geom_quantile_box(aes(fill=GROUP)) + 
    geom_boxplot(aes(col=METFORMIN), fill=NA, outlier.shape = NA) +
    geom_jitter(aes(col=METFORMIN), alpha=0.5, 
                position = position_jitterdodge()) +
    scale_fill_manual(values=c('#70737250', '#E95F0A')) + 
    scale_colour_manual(values=rev(c('#94CCFB', '#052696'))) +
    theme_publication(panel.grid = 'major_y') + 
    facet_grid(~species)
ggsave(g, filename = './figures/misc/t2d_confounding.pdf',
       width = 6, height = 4, useDingbats=FALSE)
l <- as.numeric(as.factor(meta$GROUP))
l <- ((l-2) * 2 ) + 1
names(l) <- rownames(meta)

c <- as.numeric(as.factor(meta$METFORMIN))
names(c) <- rownames(meta)
c <- data.frame(conf=c)


feat.x <- feat[x,names(l)]


p.limma <-test.via.limma(data=asin(sqrt(feat.x)), label=l, conf=NULL)
p.lm <-test.lm(data=log10(feat.x+1e-05), label=l, conf=NULL)
p.wilcox <- test.wilcoxon(data=feat.x, label=l, conf=NULL)

p.limma.c <-test.via.limma(data=asin(sqrt(feat.x)), label=l, conf=c)
p.lm.c <-test.lme(data=log10(feat.x+1e-05), label=l, conf=c)
p.wilcox.c <- test.wilcoxon(data=feat.x, label=l, conf=c)

enframe(p.limma) %>% 
  mutate(test='limma') %>% 
  mutate(corrected=FALSE) %>% 
  bind_rows(enframe(p.limma.c) %>% 
              mutate(test='limma') %>% 
              mutate(corrected=TRUE)) %>% 
  bind_rows(enframe(p.lm) %>% 
              mutate(test='lm') %>% 
              mutate(corrected=FALSE)) %>% 
  bind_rows(enframe(p.lm.c) %>% 
              mutate(test='lm') %>% 
              mutate(corrected=TRUE)) %>% 
  bind_rows(enframe(p.wilcox) %>% 
              mutate(test='wilcoxon') %>% 
              mutate(corrected=FALSE)) %>% 
  bind_rows(enframe(p.wilcox.c) %>% 
              mutate(test='wilcoxon') %>% 
              mutate(corrected=TRUE)) %>% 
  mutate(label=sprintf(fmt='%.1E', value)) %>%
  ggplot(aes(x=test, y=corrected, fill=test)) + 
  geom_tile(aes(fill=test)) + 
  geom_text(aes(label=label), colour='white') +
  scale_fill_manual(values=unlist(test.colours)) +
  facet_grid(~name)
    
