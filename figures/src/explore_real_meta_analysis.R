# ##############################################################################
#
## Look at confounder-stuff in the CRC/CD data
#
# ##############################################################################

library("tidyverse")
library("here")
library("coin")
library("ggembl")
library("cowplot")
library("progress")
library('lmerTest')
# ##############################################################################
# function

# gFC plot
.f_contrast_gfc <- function(tbl.feat, tbl.meta, label, confounder){
  stopifnot(label %in% colnames(tbl.meta))
  stopifnot(confounder %in% colnames(tbl.meta))
  stopifnot("Sample_ID" %in% colnames(tbl.meta))
  stopifnot(all(tbl.meta$Sample_ID %in% colnames(tbl.feat)))
  
  stopifnot(length(na.omit(unique(tbl.meta[[label]])))==2)
  stopifnot(length(na.omit(unique(tbl.meta[[confounder]])))==2)
  l.a <- na.omit(unique(tbl.meta[[label]]))[1]
  l.b <- na.omit(unique(tbl.meta[[label]]))[2]
  c.a <- na.omit(unique(tbl.meta[[confounder]]))[1]
  c.b <- na.omit(unique(tbl.meta[[confounder]]))[2]
  
  # prevalence filter
  tbl.feat <- tbl.feat[,tbl.meta$Sample_ID]
  tbl.feat <- tbl.feat[rowMeans(tbl.feat!=0) > 0.05,]
  
  # compute gFC with the label
  gfc.label <- vapply(rownames(tbl.feat), FUN = function(x){
    q.a <- quantile(log10(tbl.feat[x,tbl.meta %>% 
                                     filter(!!as.symbol(label)==l.a) %>%
                                     pull(Sample_ID)] + 1e-05), 
                    probs=seq(0.1, 0.9, by=0.05))
    q.b <- quantile(log10(tbl.feat[x,tbl.meta %>% 
                                     filter(!!as.symbol(label)==l.b) %>%
                                     pull(Sample_ID)] + 1e-05), 
                    probs=seq(0.1, 0.9, by=0.05))
    mean(q.a-q.b)
  }, FUN.VALUE = double(1))
  
  # compute gFC with the confounder
  gfc.confounder <- vapply(rownames(tbl.feat), FUN = function(x){
    q.a <- quantile(log10(tbl.feat[x,tbl.meta %>% 
                                     filter(!!as.symbol(confounder)==c.a) %>%
                                     pull(Sample_ID)] + 1e-05), 
                    probs=seq(0.1, 0.9, by=0.05))
    q.b <- quantile(log10(tbl.feat[x,tbl.meta %>% 
                                     filter(!!as.symbol(confounder)==c.b) %>%
                                     pull(Sample_ID)] + 1e-05), 
                    probs=seq(0.1, 0.9, by=0.05))
    mean(q.a-q.b)
  }, FUN.VALUE = double(1))
  
  enframe(gfc.label, value='label') %>% 
    full_join(enframe(gfc.confounder, value='confounder'), by='name')
}

# test normal/blocked
.f_test_naive_blocked <- function(tbl.feat, tbl.meta, label, confounder,
                                  type='wilcoxon'){
  stopifnot(label %in% colnames(tbl.meta))
  stopifnot(confounder %in% colnames(tbl.meta))
  stopifnot("Sample_ID" %in% colnames(tbl.meta))
  stopifnot(all(tbl.meta$Sample_ID %in% colnames(tbl.feat)))
  
  stopifnot(length(na.omit(unique(tbl.meta[[label]])))==2)
  
  # prevalence filter
  tbl.feat <- tbl.feat[,tbl.meta$Sample_ID]
  tbl.feat <- tbl.feat[rowMeans(tbl.feat!=0) > 0.05,]
  
  # compute gFC with the label
  if (type=='wilcoxon'){
    p.val.naive <- vapply(rownames(tbl.feat), FUN = function(x){
      t <- wilcox.test(tbl.feat[x,tbl.meta$Sample_ID]~tbl.meta[[label]])
      t$p.value
    }, FUN.VALUE = double(1))
  
    pb <- progress_bar$new(total=nrow(tbl.feat))
    p.val.blocked <- vapply(rownames(tbl.feat), FUN = function(x){
      t <- wilcox_test(tbl.feat[x,tbl.meta$Sample_ID]~
                         as.factor(tbl.meta[[label]])|
                         as.factor(tbl.meta[[confounder]]))
      pb$tick()
      pvalue(t)
    }, FUN.VALUE = double(1))
  } else if (type=='lme'){
    p.val.naive <- vapply(rownames(tbl.feat), FUN = function(x){
      fit <- lm(log10(tbl.feat[x,tbl.meta$Sample_ID] + 1e-05)~
                  tbl.meta[[label]])
      res <- coefficients(summary(fit))
      res[2,4]
    }, FUN.VALUE = double(1))
    
    pb <- progress_bar$new(total=nrow(tbl.feat))
    p.val.blocked <- vapply(rownames(tbl.feat), FUN = function(x){
      df <- tbl.meta
      df[['taxon']] <- log10(tbl.feat[x,tbl.meta$Sample_ID] + 1e-05)
      fit <- suppressMessages(
        lmer(data=df,  paste0('taxon~', label, '+(1|', confounder, ')')))
      res <- coefficients(summary(fit))
      pb$tick()
      res[2,5]
    }, FUN.VALUE = double(1))
  } else if (type=='limma') {
    x <- as.data.frame(tbl.meta)
    x.red <- x[,which(colnames(x) %in% c(label, confounder))]
    x.red <- vapply(colnames(x.red), 
                    FUN = function(x){as.numeric(as.factor(x.red[[x]]))}, 
                    FUN.VALUE = double(nrow(x.red)))
    rownames(x.red) <- x$Sample_ID
    fit <- limma::lmFit(tbl.feat, 
                        design = x.red[,which(colnames(x.red)==label),
                                       drop=FALSE])
    res <- limma::eBayes(fit)
    p.val.naive <- res$p.value[,label]
    
    fit.blocked <- limma::lmFit(tbl.feat, design = x.red)
    res.blocked <- limma::eBayes(fit.blocked)
    p.val.blocked <- res.blocked$p.value[,label]
  } else {
    stop("unknown test type")
  }
  
  enframe(p.val.naive, value='naive') %>% 
    full_join(enframe(p.val.blocked, value='blocked'), by='name') %>% 
    mutate(adj.naive=p.adjust(naive, method='fdr')) %>% 
    mutate(adj.blocked=p.adjust(blocked, method='fdr')) %>% 
    mutate(found.naive=adj.naive < 0.05, found.blocked=adj.blocked < 0.05)
  
}

# ##############################################################################
# load data
meta.crc <- read_tsv(here('data', 'meta_crc.tsv'))
feat.crc <- read.table(here('data', 'motus_crc_meta.tsv'), sep='\t',
                       stringsAsFactors = FALSE, check.names = FALSE, 
                       quote = '', comment.char = '')
feat.crc <- prop.table(as.matrix(feat.crc), 2)

# plot study-confounding numbers  
g.study.crc <- meta.crc %>% 
  ggplot(aes(x=Study, fill=Group)) + 
    geom_bar() + 
    theme_publication() + 
    scale_fill_manual(values=c('#18974C', '#A8A99E'))

# same for CD
meta.cd <- read_tsv('../../SIAMCAT/siamcat_paper/ibd_meta_analysis/data/meta_all.tsv')
meta.cd.ind <- meta.cd %>% 
  group_by(Individual_ID) %>% 
  filter(Timepoint==min(Timepoint)) %>% 
  ungroup()
feat.cd <- read.table('../../SIAMCAT/siamcat_paper/ibd_meta_analysis/data/feat_genus.tsv', 
                      sep='\t', stringsAsFactors = FALSE, check.names = FALSE)
feat.cd <- as.matrix(feat.cd)

g.study.cd <- meta.cd.ind %>% 
  ggplot(aes(x=Study, fill=Group)) + 
  geom_bar() + 
  theme_publication() + 
  scale_fill_manual(values=c('#18974C', '#A8A99E'))


# ##############################################################################
# loop over all comparisons
.f_compare_tests <- function(res.now, res.all){
  x <- res.now %>% 
    select(name, found.naive, found.blocked) %>% 
    full_join(res.all %>% transmute(name, all=found.blocked), by='name')
  x.naive <- table(x$found.naive, x$all)
  x.blocked <- table(x$found.blocked, x$all)
  
  fpr.naive <- x.naive[2,1]/sum(x.naive[,1])
  fpr.blocked <- x.blocked[2,1]/sum(x.blocked[,1])
  
  miss.rate.naive <- x.naive[1,2]/sum(x.naive[,2])
  miss.rate.blocked <- x.blocked[1,2]/sum(x.blocked[,2])
  
  return(list('fpr.blocked'=fpr.blocked, 'fpr.naive'=fpr.naive,
              'miss.rate.blocked'=miss.rate.blocked, 
              'miss.rate.naive'=miss.rate.naive))
}

fn.fpr <- here('figures','figure_conf_study', 'false_predictions.tsv')
if (!file.exists(fn.fpr)){
  df.plot.v <- tibble(comparison=character(0), disease=character(0), 
                      n.samples=double(0), test=character(0),
                      miss.rate=double(0), FPR=double(0), v.value=double(0))
  studies.crc <- unique(meta.crc$Study)
  studies.cd <- unique(meta.cd.ind$Study)
  study.selection <- expand_grid('a'=c(0,1), 'b'=c(0,1), 'v'=c(0,1), 
                                 'd'=c(0,1), 'e'=c(0,1))
  # CRC all studies
  crc.all.w <- .f_test_naive_blocked(feat.crc, meta.crc, 
                                   'Group', 'Study', 'wilcoxon')
  crc.all.lme <- .f_test_naive_blocked(feat.crc, meta.crc, 
                                       'Group', 'Study', 'lme') 
  crc.all.li <- .f_test_naive_blocked(feat.crc, meta.crc, 
                                       'Group', 'Study', 'limma') 
  # CD all studies
  cd.all.w <- .f_test_naive_blocked(feat.cd, meta.cd.ind, 
                                     'Group', 'Study', 'wilcoxon') 
  cd.all.lme <- .f_test_naive_blocked(feat.cd, meta.cd.ind, 
                                       'Group', 'Study', 'lme') 
  cd.all.li <- .f_test_naive_blocked(feat.cd, meta.cd.ind, 
                                      'Group', 'Study', 'limma') 
  for (i in seq_len(nrow(study.selection))){
    idx <- unlist(study.selection[i,])
    if (sum(idx) < 2) next()
    
    # first CRC
    message(i, '-CRC')
    sel.studies <- studies.crc[idx==1]
    
    meta.red <- meta.crc %>% 
      filter(Study %in% sel.studies)
    # compute cramers v
    phi <- lsr::cramersV(meta.red$Group, meta.red$Study)
    # wilcoxon test & LME & limma
    df.test.w <- .f_test_naive_blocked(feat.crc, meta.red, 
                                       'Group', 'Study', 'wilcoxon') 
    w <- .f_compare_tests(df.test.w, crc.all.w)
    df.plot.v <- df.plot.v %>% 
      add_row(comparison=paste(sel.studies, collapse = ','),
              disease='CRC', n.samples=nrow(meta.red),
              miss.rate=c(w$miss.rate.blocked, w$miss.rate.naive), 
              FPR=c(w$fpr.blocked, w$fpr.naive), 
              v.value=phi, test=c('wilcoxon_blocked', 'wilcoxon'))
    # LME
    df.test.lme <- .f_test_naive_blocked(feat.crc, meta.red, 
                                       'Group', 'Study', 'lme') 
    lme <- .f_compare_tests(df.test.lme, crc.all.lme)
    df.plot.v <- df.plot.v %>% 
      add_row(comparison=paste(sel.studies, collapse = ','),
              disease='CRC', n.samples=nrow(meta.red),
              miss.rate=c(lme$miss.rate.blocked, lme$miss.rate.naive), 
              FPR=c(lme$fpr.blocked, lme$fpr.naive), 
              v.value=phi, test=c('LM_blocked', 'LM'))
    # limma
    df.test.limma <- .f_test_naive_blocked(feat.crc, meta.red, 
                                         'Group', 'Study', 'limma') 
    limma <- .f_compare_tests(df.test.limma, crc.all.li)
    df.plot.v <- df.plot.v %>% 
      add_row(comparison=paste(sel.studies, collapse = ','),
              disease='CRC', n.samples=nrow(meta.red),
              miss.rate=c(limma$miss.rate.blocked, limma$miss.rate.naive), 
              FPR=c(limma$fpr.blocked, limma$fpr.naive), 
              v.value=phi, test=c('limma_blocked', 'limma'))
    # second CD
    message(i, '-CD')
    sel.studies <- studies.cd[idx==1]
    
    meta.red <- meta.cd.ind %>% 
      filter(Study %in% sel.studies)
    # compute cramers v
    phi <- lsr::cramersV(meta.red$Group, meta.red$Study)
    # wilcoxon test & LME & limma
    df.test.w <- .f_test_naive_blocked(feat.cd, meta.red, 
                                       'Group', 'Study', 'wilcoxon') 
    w <- .f_compare_tests(df.test.w, cd.all.w)
    df.plot.v <- df.plot.v %>% 
      add_row(comparison=paste(sel.studies, collapse = ','),
              disease='CD', n.samples=nrow(meta.red),
              miss.rate=c(w$miss.rate.blocked, w$miss.rate.naive), 
              FPR=c(w$fpr.blocked, w$fpr.naive), 
              v.value=phi, test=c('wilcoxon_blocked', 'wilcoxon'))
    # LME
    df.test.lme <- .f_test_naive_blocked(feat.cd, meta.red, 
                                         'Group', 'Study', 'lme') 
    lme <- .f_compare_tests(df.test.lme, cd.all.lme)
    df.plot.v <- df.plot.v %>% 
      add_row(comparison=paste(sel.studies, collapse = ','),
              disease='CD', n.samples=nrow(meta.red),
              miss.rate=c(lme$miss.rate.blocked, lme$miss.rate.naive), 
              FPR=c(lme$fpr.blocked, lme$fpr.naive), 
              v.value=phi, test=c('LM_blocked', 'LM'))
    # limma
    df.test.limma <- .f_test_naive_blocked(feat.cd, meta.red, 
                                           'Group', 'Study', 'limma') 
    limma <- .f_compare_tests(df.test.limma, cd.all.li)
    df.plot.v <- df.plot.v %>% 
      add_row(comparison=paste(sel.studies, collapse = ','),
              disease='CD', n.samples=nrow(meta.red),
              miss.rate=c(limma$miss.rate.blocked, limma$miss.rate.naive), 
              FPR=c(limma$fpr.blocked, limma$fpr.naive), 
              v.value=phi, test=c('limma_blocked', 'limma'))
    
  }
  write_tsv(df.plot.v, file = fn.fpr)
} else {
  df.plot.v <- read_tsv(fn.fpr)
}
g <- df.plot.v %>% 
  pivot_longer(cols=c(miss.rate, FPR)) %>% 
  mutate(blocked=str_detect(test, 'blocked')) %>% 
  mutate(test=str_remove(test, '_blocked')) %>% 
  filter(name=='FPR') %>% 
  mutate(label=case_when(v.value > 0~comparison, TRUE~'')) %>% 
  ggplot(aes(x=v.value, y=value, col=blocked)) + 
    geom_point(aes(shape=disease)) + 
    facet_grid(~test) +
    theme_publication(panel.grid='major') + 
    geom_smooth(method='lm') + 
    xlab("Cramers's V") + 
    ylab('False positive rate (compared to all data combined)') + 
    scale_colour_manual(values=c('#D41645', '#3B6FB6')) + 
    ggrepel::geom_text_repel(aes(label=label), size=1)
    
ggsave(g, filename = here('figures', 'figure_conf_study', 'cramers.pdf'),
       width = 7, height = 4)


# 





# comparison of study/group effects in pairwise comparisons
pdf(here('figures', 'figure_conf_study', 'real_comparisons_CRC.pdf'),
    width = 5, height = 3, useDingbats = FALSE)
studies <- unique(meta.crc$Study)
df.plot.phi <- tibble(Study1=character(0), Study2=character(0),
                      miss.rate=double(0), FPR=double(0), phi=double(0))
for (i in seq_along(studies)){
  s1 <- studies[i]
  if (i==length(studies)) next()
  for (i2 in (i+1):length(studies)){
    s2 <- studies[i2]
    message(s1, '-', s2)
    meta.red <- meta.crc %>%  filter(Study%in%c(s1, s2))
    df.plot <- .f_contrast_gfc(feat.crc, meta.red,'Group', 'Study')
    df.test.w <- .f_test_naive_blocked(feat.crc, meta.red, 'Group', 'Study', 'wilcoxon')
    df.test.l <- .f_test_naive_blocked(feat.crc, meta.red, 'Group', 'Study', 'lme')
    g1 <- df.plot %>% 
      ggplot(aes(x=abs(label), y=abs(confounder))) + 
      geom_abline(slope = 1, intercept = 0) + 
      geom_point() + 
      xlim(0,2) + ylim(0,2) + 
      theme_publication() + 
      ggtitle(paste0(s1, '-', s2))
    sp <- cor(df.test$blocked, df.test$naive, method='spearman')
    df.test <- df.test %>% 
      mutate(adj.naive=p.adjust(naive, method='fdr')) %>% 
      mutate(adj.blocked=p.adjust(blocked, method='fdr')) %>% 
      mutate(found.naive=adj.naive < 0.05, found.blocked=adj.blocked < 0.05)
    
    # sensitivity/specificity?
    x <- table(df.test$found.naive, df.test$found.blocked)
    miss.rate <- x[1,2]/sum(x[,2])
    fpr <- x[2,1]/sum(x[,1])
    
    g2 <- df.test %>% 
      ggplot(aes(x=-log10(naive), y=-log10(blocked))) +
      geom_abline(slope = 1, intercept = 0) + 
      geom_point() + 
      theme_publication() +
      annotate(geom='text', x=5, y=10, 
               label=paste0('rho=', sprintf(fmt='%.2f', sp)))
    g <- plot_grid(g1, g2, nrow = 1)
    print(g)
    
    df.plot.phi <- df.plot.phi %>% 
      add_row(Study1=s1, Study2=s2, miss.rate=miss.rate, 
              FPR=fpr, phi=-phi(table(meta.red$Group, meta.red$Study)))
  }
}
dev.off()

# ##############################################################################
# other confounders? diet/vegetarian_diet/sampling_relative_to_colonnoscopy

table(meta.crc$Study, meta.crc$Vegetarian)
table(meta.crc$Study, meta.crc$Diabetes)
table(meta.crc$Study, meta.crc$Sampling_rel_to_colonoscopy)


df.plot <- .f_contrast_gfc(feat.crc, 
                           meta.crc %>% filter(Study=='AT-CRC'),
                           'Group', 'Vegetarian')
df.test.v <- .f_test_naive_blocked(feat.crc, 
                                 meta.crc %>% filter(Study=='AT-CRC'),
                                 'Vegetarian', 'Group')

df.plot <- .f_contrast_gfc(feat.crc, 
                           meta.crc %>% filter(Study=='AT-CRC'),
                           'Group', 'Diabetes')
df.test <- .f_test_naive_blocked(feat.crc, 
                                 meta.crc %>% filter(Study=='AT-CRC'),
                                 'Group', 'Diabetes')

df.plot <- .f_contrast_gfc(feat.crc, 
                           meta.crc %>% filter(Study=='CN-CRC'),
                           'Group', 'Diabetes')
df.test <- .f_test_naive_blocked(feat.crc, 
                                 meta.crc %>% filter(Study=='CN-CRC'),
                                 'Group', 'Diabetes')

df.plot <- .f_contrast_gfc(feat.crc, 
                           meta.crc %>% filter(Study=='DE-CRC'),
                           'Group', 'Diabetes')
df.test <- .f_test_naive_blocked(feat.crc, 
                                 meta.crc %>% filter(Study=='DE-CRC'),
                                 'Group', 'Diabetes')

df.plot <- .f_contrast_gfc(feat.crc, 
                           meta.crc %>% filter(Study=='CN-CRC'),
                           'Group', 'Sampling_rel_to_colonoscopy')
df.test <- .f_test_naive_blocked(feat.crc, 
                                 meta.crc %>% filter(Study=='CN-CRC'),
                                 'Group', 'Sampling_rel_to_colonoscopy')

df.plot %>% 
  ggplot(aes(x=abs(label), y=abs(confounder))) + 
    geom_abline(slope = 1, intercept = 0) + 
    geom_point() + 
    xlim(0,2) + ylim(0,2) + 
    theme_publication()

df.test %>% 
  ggplot(aes(x=-log10(naive), y=-log10(blocked))) +
    geom_abline(slope = 1, intercept = 0) + 
    geom_point() + 
    theme_publication()

# artificial bias
meta.crc.temp <- bind_rows(
  meta.crc %>% 
    filter(Study=='CN-CRC') %>% filter(Group=='CRC') %>% slice_sample(n=40),
  meta.crc %>% 
    filter(Study=='CN-CRC') %>% filter(Group=='CTR') %>% slice_sample(n=20),
  meta.crc %>% 
    filter(Study=='AT-CRC') %>% filter(Group=='CRC') %>% slice_sample(n=20),
  meta.crc %>% 
    filter(Study=='AT-CRC') %>% filter(Group=='CTR') %>% slice_sample(n=40))

df.plot <- .f_contrast_gfc(feat.crc, meta.crc.temp, 'Group', 'Study')
df.test <- .f_test_naive_blocked(feat.crc, meta.crc.temp, 'Group', 'Study')

df.test.all <- .f_test_naive_blocked(
  feat.crc, meta.crc %>% 
    filter(!is.na(Sampling_rel_to_colonoscopy)) %>% 
    mutate(conf=case_when(Study!='CN-CRC'~Study, 
                          TRUE~paste0(Study, '-', 
                                      Sampling_rel_to_colonoscopy))),
  'Group', 'conf') %>% 
  mutate(p.adj.naive=p.adjust(naive, method='fdr'),
         p.adj.block=p.adjust(blocked, method='fdr'))
df.test.all %>% 
  ggplot(aes(x=-log10(p.adj.naive), y=-log10(p.adj.block))) +
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() + 
  theme_publication()
  
# ##############################################################################
# CD-stuff?

meta.cd <- read_tsv('../../SIAMCAT/siamcat_paper/ibd_meta_analysis/data/meta_all.tsv')
meta.cd.ind <- meta.cd %>% 
  group_by(Individual_ID) %>% 
  filter(Timepoint==min(Timepoint)) %>% 
  ungroup()
feat.cd <- read.table('../../SIAMCAT/siamcat_paper/ibd_meta_analysis/data/feat_genus.tsv', 
                      sep='\t', stringsAsFactors = FALSE, check.names = FALSE)
feat.cd <- as.matrix(feat.cd)

g.study.2 <- meta.cd.ind %>% 
  ggplot(aes(x=Study, fill=Group)) + 
  geom_bar() + 
  theme_publication() + 
  scale_fill_manual(values=c('#18974C', '#A8A99E'))


df.test.all <- .f_test_naive_blocked(feat.cd, meta.cd.ind, 'Group', 'Study') %>% 
  mutate(p.adj.naive=p.adjust(naive, method='fdr'),
         p.adj.block=p.adjust(blocked, method='fdr'))

pdf(here('figures', 'figure_conf_study', 'real_comparisons_CD.pdf'),
    width = 5, height = 3, useDingbats = FALSE)
studies <- unique(meta.cd.ind$Study)
for (i in seq_along(studies)){
  s1 <- studies[i]
  if (i==length(studies)) next()
  for (i2 in (i+1):length(studies)){
    s2 <- studies[i2]
    meta.red <- meta.cd.ind %>% 
      filter(Study%in%c(s1, s2))
    df.plot <- .f_contrast_gfc(feat.cd, meta.red, 'Group', 'Study')
    df.test <- .f_test_naive_blocked(feat.cd, meta.red, 'Group', 'Study')
    g1 <- df.plot %>% 
      ggplot(aes(x=abs(label), y=abs(confounder))) + 
      geom_abline(slope = 1, intercept = 0) + 
      geom_point() + 
      xlim(0,2) + ylim(0,2) + 
      theme_publication()
    sp <- cor(df.test$blocked, df.test$naive, method='spearman')
    df.test <- df.test %>% 
      mutate(adj.naive=p.adjust(naive, method='fdr')) %>% 
      mutate(adj.blocked=p.adjust(blocked, method='fdr')) %>% 
      mutate(found.naive=adj.naive < 0.05, found.blocked=adj.blocked < 0.05)
    
    # sensitivity/specificity?
    x <- table(df.test$found.naive, df.test$found.blocked)
    miss.rate <- x[1,2]/sum(x[,2])
    fpr <- x[2,1]/sum(x[,1])
    g2 <- df.test %>% 
      ggplot(aes(x=-log10(naive), y=-log10(blocked))) +
      geom_abline(slope = 1, intercept = 0) + 
      geom_point() + 
      theme_publication() + 
      annotate(geom='text', x=1, y=10, 
               label=paste0('rho=', sprintf(fmt='%.2f', sp)))
    g <- plot_grid(g1, g2, nrow = 1)
    print(g)
    df.plot.phi <- df.plot.phi %>% 
      add_row(Study1=s1, Study2=s2, miss.rate=miss.rate, 
              FPR=fpr, phi=-phi(table(meta.red$Group, meta.red$Study)))
  }
}
dev.off()



df.plot %>% mutate(l.diff=(abs(label))) %>% 
  select(name, l.diff) %>% full_join(
  df.test %>% mutate(naive=-log10(naive), blocked=-log10(blocked)) %>% 
    mutate(diff=blocked-naive) %>% select(name, diff), by='name') %>% 
  ggplot(aes(x=l.diff, y=diff)) + 
    geom_point()

df.test.all <- .f_test_naive_blocked(
  feat.cd, meta.cd, 'Group', 'Study') %>% 
  mutate(p.adj.naive=p.adjust(naive, method='fdr'),
         p.adj.block=p.adjust(blocked, method='fdr'))
df.test.all %>% 
  ggplot(aes(x=-log10(p.adj.naive), y=-log10(p.adj.block))) +
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() + 
  theme_publication()


# save
g.study <- cowplot::plot_grid(g.study.1, g.study.2, ncol = 1)
ggsave(g.study, filename = here('figures','figure_conf_study',
                                'real_study_confounding.pdf'),
       width = 5, height = 8, useDingbats=FALSE)

# some examples

g.exp.1 <- .f_contrast_gfc(feat.crc, 
                           meta.crc %>% 
                             filter(Study%in%c('CN-CRC', 'AT-CRC')), 
                           'Group', 'Study') %>% 
  ggplot(aes(x=abs(label), y=abs(confounder))) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() + 
  xlim(0,2) + ylim(0,2) + 
  theme_publication() + 
  ggtitle('CN-vs-AT')

g.exp.2 <- .f_contrast_gfc(feat.crc, 
                           meta.crc %>% 
                             filter(Study%in%c('AT-CRC', 'FR-CRC')), 
                           'Group', 'Study') %>% 
  ggplot(aes(x=abs(label), y=abs(confounder))) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() + 
  xlim(0,2) + ylim(0,2) + 
  theme_publication() + 
  ggtitle('AT-vs-FR')

g.exp.3 <- .f_contrast_gfc(feat.cd, 
                           meta.cd.ind %>% 
                             filter(Study%in%c('Lewis_2015', 'metaHIT')), 
                           'Group', 'Study') %>% 
  ggplot(aes(x=abs(label), y=abs(confounder))) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() + 
  xlim(0,2) + ylim(0,2) + 
  theme_publication() + 
  ggtitle('lewis_Metahit')

g.exp.4 <- .f_contrast_gfc(feat.cd, 
                           meta.cd.ind %>% 
                             filter(Study%in%c('Franzosa_2019', 'He_2017')), 
                           'Group', 'Study') %>% 
  ggplot(aes(x=abs(label), y=abs(confounder))) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() + 
  xlim(0,2) + ylim(0,2) + 
  theme_publication() + 
  ggtitle('HE-Franzosa')

g.all <- cowplot::plot_grid(g.exp.1, g.exp.2, g.exp.3, g.exp.4, nrow = 2)
ggsave(g.all, filename = here('figures','figure_conf_study',
                                'real_study_confounding_gFC.pdf'),
       width = 8, height = 8, useDingbats=FALSE)




# ##############################################################################
# ICI stuff?




# load data
data.location <- 'https://intranet.embl.de/download/zeller/'

# metadata
meta.routy <- read_tsv(paste0(data.location, 'metadata/meta_Routy.tsv'))
meta.matson <- read_tsv(paste0(data.location, 'metadata/meta_Matson.tsv'))
# motus
feat.ab.routy <- read.table(paste0(data.location, 
                                   'tax_profiles/routy_motus.v2.1.motus'),
                            stringsAsFactors = FALSE, check.names = FALSE,
                            quote = '', comment.char = '#', sep='\t',
                            header = TRUE, row.names = 1)
feat.rel.routy <- prop.table(as.matrix(feat.ab.routy), 2)
feat.ab.matson <- read.table(paste0(data.location, 
                                    'tax_profiles/matson_motus.v2.1.motus'),
                             stringsAsFactors = FALSE, check.names = FALSE,
                             quote = '', comment.char = '#', sep='\t',
                             header = TRUE, row.names = 1)
feat.rel.matson <- prop.table(as.matrix(feat.ab.matson), 2)

feat.ici <- cbind(feat.rel.routy, feat.rel.matson)
meta.ici <- meta.matson %>% 
  mutate(Study='Matson') %>% 
  mutate(Disease='MEL') %>% 
  bind_rows(meta.routy %>% 
              group_by(Individual_ID) %>% 
              filter(Timepoint==min(Timepoint)) %>% 
              select(Sample_ID, Individual_ID, Group, Disease) %>% 
              mutate(Study='Routy') %>% 
              mutate(Group=case_when(Group=='Dead'~'no response',
                                     Group=='Progression'~'no response',
                                     TRUE~'response')))
df.plot <- .f_contrast_gfc(feat.ici, meta.ici, 'Group', 'Study')
df.test <- .f_test_naive_blocked(feat.ici, meta.ici, 'Group', 'Study')
df.plot %>% 
  ggplot(aes(x=label, y=confounder)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() + 
  xlim(-2,2) + ylim(-2,2) + 
  theme_publication()
df.test %>% 
  full_join(df.plot, by='name') %>% 
  mutate(naive=-log10(naive), blocked=-log10(blocked)) %>% 
  mutate(naive=naive*sign(label), blocked=blocked*sign(label)) %>% 
  ggplot(aes(x=naive, y=blocked)) +
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() + 
  theme_publication()  
