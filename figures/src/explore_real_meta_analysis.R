# ##############################################################################
#
## Look at real meta-analysis setting for CRC/CD data
##
##  This involves looking at Study as a confounder
##  Also includes testing all association models on the combined data
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
  df.temp <- t(tbl.feat)
  df.temp <- as.data.frame(df.temp[tbl.meta$Sample_ID,])
  df.temp$group <- as.numeric(as.factor(tbl.meta[[label]]))
  # compute gFC with the label
  if (type=='wilcoxon'){
    p.val.naive <- vapply(rownames(tbl.feat), FUN = function(x){
      t <- wilcox.test(df.temp[[x]]~df.temp[["group"]])
      t$p.value
    }, FUN.VALUE = double(1))
  
    pb <- progress_bar$new(total=nrow(tbl.feat))
    p.val.blocked <- vapply(rownames(tbl.feat), FUN = function(x){
      t <- wilcox_test(df.temp[[x]]~
                         as.factor(tbl.meta[[label]])|
                         as.factor(tbl.meta[[confounder]]))
      pb$tick()
      pvalue(t)
    }, FUN.VALUE = double(1))
  } else if (type=='lme'){
    p.val.naive <- vapply(rownames(tbl.feat), FUN = function(x){
      fit <- lm(log10(df.temp[[x]] + 1e-05)~
                  tbl.meta[[label]])
      res <- coefficients(summary(fit))
      res[2,4]
    }, FUN.VALUE = double(1))
    
    pb <- progress_bar$new(total=nrow(tbl.feat))
    p.val.blocked <- vapply(rownames(tbl.feat), FUN = function(x){
      df <- tbl.meta
      df[['taxon']] <- log10(df.temp[[x]]+ 1e-05)
      if (length(unique(df$taxon)) == 1){
        pb$tick()
        return(1)
      }
      fit <- suppressMessages(
        lmer(data=df,  paste0('taxon~', label, '+(1|', confounder, ')')))
      res <-  coefficients(summary(fit))
      pb$tick()
      return(res[2,5])
    }, FUN.VALUE = double(1))
  } else if (type=='limma') {
    # browser()
    tbl.feat <- log10(tbl.feat + 1e-05)
    x <- as.data.frame(tbl.meta)
    x.red <- x[,which(colnames(x) %in% c(label, confounder))]
    x.red <- vapply(colnames(x.red), 
                    FUN = function(x){as.numeric(as.factor(x.red[[x]]))}, 
                    FUN.VALUE = double(nrow(x.red)))
    rownames(x.red) <- x$Sample_ID
    design = data.frame(
      intercept=-1, 
      label=x.red[,which(colnames(x.red)==label),
                  drop=FALSE])
    fit <- limma::lmFit(tbl.feat, design = design)
    res <- limma::eBayes(fit)
    p.val.naive <- res$p.value[,label]
    
    fit.blocked <- limma::lmFit(tbl.feat, design = design, 
                                block = x.red[,confounder], correlation=0.5)
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
test.colours <- yaml::read_yaml(here('files', 'test_colours.yml'))

# ##############################################################################
# load data
meta.crc <- read_tsv(here('data', 'meta_crc.tsv'))
feat.crc <- read.table(here('data', 'motus_crc_rel_meta.tsv'), sep='\t',
                       stringsAsFactors = FALSE, check.names = FALSE, 
                       quote = '', comment.char = '')
feat.crc <- as.matrix(feat.crc)

# plot study-confounding numbers  
g.study.crc <- meta.crc %>% 
  ggplot(aes(x=Study, fill=Group)) + 
    geom_bar() + 
    theme_publication() + 
    scale_fill_manual(values=c('#18974C', '#A8A99E'))
g.study.crc.rel <- meta.crc %>% 
  group_by(Study, Group) %>% 
  tally() %>% 
  group_by(Study) %>% 
  mutate(prop=n/sum(n)) %>% 
  ggplot(aes(x=Study, fill=Group, y=prop)) + 
    geom_col() + 
    theme_publication() + 
    scale_fill_manual(values=c('#18974C', '#A8A99E')) + 
    ylab("Proportion")
ggsave(g.study.crc, 
       filename = here('figures', 'figure_conf_study', 'CRC_studies.pdf'),
       width = 4, height = 4)
ggsave(g.study.crc.rel, 
       filename = here('figures', 'figure_conf_study', 'CRC_studies_rel.pdf'),
       width = 4, height = 4)

# same for CD
meta.cd <- read_tsv(here('data', 'meta_cd.tsv'))
feat.cd <- read.table(here('data', 'motus_cd_rel_meta.tsv'), sep='\t',
                      stringsAsFactors = FALSE, check.names = FALSE, 
                      quote = '', comment.char = '')
feat.cd <- as.matrix(feat.cd)

g.study.cd <- meta.cd %>% 
  ggplot(aes(x=Study, fill=Group)) + 
  geom_bar() + 
  theme_publication() + 
  scale_fill_manual(values=c('#18974C', '#A8A99E'))
g.study.cd.rel <- meta.cd %>% 
  group_by(Study, Group) %>% 
  tally() %>% 
  group_by(Study) %>% 
  mutate(prop=n/sum(n)) %>% 
  ggplot(aes(x=Study, fill=Group, y=prop)) + 
    geom_col() + 
    theme_publication() + 
    scale_fill_manual(values=c('#18974C', '#A8A99E')) + 
    ylab("Proportion")
ggsave(g.study.cd, 
       filename = here('figures', 'figure_batch_effects', 'CD_studies.pdf'),
       width = 4, height = 4)
ggsave(g.study.cd.rel, 
       filename = here('figures', 'figure_batch_effects', 'CD_studies_rel.pdf'),
       width = 4, height = 4)

# ##############################################################################
# loop over all comparisons
.f_compare_tests <- function(res.now, res.all){
  x <- res.now %>% 
    select(name, found.naive, found.blocked) %>% 
    full_join(res.all %>% transmute(name, all=found.blocked), by='name')
  x.naive <- table(x$found.naive, x$all)
  x.blocked <- table(x$found.blocked, x$all)
  
  # add false/true row if needed!
  if (!all(colnames(x.naive) == c('FALSE', 'TRUE'))){
    x <- setdiff(c('FALSE', 'TRUE'), colnames(x.naive))
    x.naive <- cbind(x.naive, c(0,0))
    colnames(x.naive)[2] <- x
  }
  if (!all(colnames(x.blocked) == c('FALSE', 'TRUE'))){
    x <- setdiff(c('FALSE', 'TRUE'), colnames(x.blocked))
    x.blocked <- cbind(x.blocked, c(0,0))
    colnames(x.blocked)[2] <- x
  }
  if (!all(rownames(x.naive) == c('FALSE', 'TRUE'))){
    x <- setdiff(c('FALSE', 'TRUE'), rownames(x.naive))
    x.naive <- rbind(x.naive, c(0,0))
    rownames(x.naive)[2] <- x
  }
  if (!all(rownames(x.blocked) == c('FALSE', 'TRUE'))){
    x <- setdiff(c('FALSE', 'TRUE'), rownames(x.blocked))
    x.blocked <- rbind(x.blocked, c(0,0))
    rownames(x.blocked)[2] <- x
  }
  
  fpr.naive <- x.naive['TRUE','FALSE']/sum(x.naive[,'FALSE'])
  fpr.blocked <- x.blocked['TRUE','FALSE']/sum(x.blocked[,'FALSE'])
  
  miss.rate.naive <- x.naive['FALSE','TRUE']/sum(x.naive[,'TRUE'])
  miss.rate.blocked <- x.blocked['FALSE','TRUE']/sum(x.blocked[,'TRUE'])
  
  return(list('fpr.blocked'=fpr.blocked, 'fpr.naive'=fpr.naive,
              'miss.rate.blocked'=miss.rate.blocked, 
              'miss.rate.naive'=miss.rate.naive))
}

fn.fpr <- here('figures','figure_batch_effects', 'false_predictions_updated.tsv')
if (!file.exists(fn.fpr)){
  df.plot.v <- tibble(comparison=character(0), disease=character(0), 
                      n.samples=double(0), test=character(0),
                      miss.rate=double(0), FPR=double(0), v.value=double(0))
  studies.crc <- unique(meta.crc$Study)
  studies.cd <- unique(meta.cd$Study)
  study.selection <- expand_grid('a'=c(0,1), 'b'=c(0,1), 'c'=c(0,1), 
                                 'd'=c(0,1), 'e'=c(0,1))
  # CRC all studies
  crc.all.w <- .f_test_naive_blocked(feat.crc, meta.crc, 
                                   'Group', 'Study', 'wilcoxon')
  crc.all.lme <- .f_test_naive_blocked(feat.crc, meta.crc, 
                                       'Group', 'Study', 'lme') 
  crc.all.li <- .f_test_naive_blocked(feat.crc, meta.crc, 
                                       'Group', 'Study', 'limma') 
  # CD all studies
  cd.all.w <- .f_test_naive_blocked(feat.cd, meta.cd, 
                                     'Group', 'Study', 'wilcoxon') 
  cd.all.lme <- .f_test_naive_blocked(feat.cd, meta.cd, 
                                       'Group', 'Study', 'lme') 
  cd.all.li <- .f_test_naive_blocked(feat.cd, meta.cd, 
                                      'Group', 'Study', 'limma') 
  cd.studies.done <- c()
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
    if (sum(idx) < 2) next()
    sel.studies <- studies.cd[idx==1]
    
    
    meta.red <- meta.cd %>% 
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
    ggrepel::geom_text_repel(aes(label=label), size=1) + 
    ylim(0, 1)
    
ggsave(g, filename = here('figures', 'figure_batch_effects', 'cramers.pdf'),
       width = 7, height = 4)


# ##############################################################################
# Look at the p-values from all association models
# to show elevated detection for some methods with poor FPR control

# if (!file.exists(here('files', 'all_tests_CD.Rdata'))){
#   stop("All tests have not yet been completed for the CD data!")
# }
if (!file.exists(here('files', 'all_tests_CRC.Rdata'))){
  stop("All tests have not yet been completed for the CRC data!")
}
# all.results.cd <- get(load(here('files', 'all_tests_CD.Rdata')))
all.results.crc <- get(load(here('files', 'all_tests_CRC.Rdata')))

# number of detection per study




# check similarity with holdout-set
studies <- meta.crc %>% 
  pull(Study) %>% unique
studies <- studies[1:5]

n <- 15
res <- list()
for (x in studies){
  
  discovery.set <- all.results.crc[[paste(setdiff(studies, x), collapse='-')]]
  validation.set <- all.results.crc[[x]]
  
  # concurrence at the top
  c10 <- map(colnames(discovery.set), .f = function(test){
    d <- sort(discovery.set[,test])
    v <- sort(validation.set[,test])
    d.bin <- names(discovery.set[,test]) %in% names(d[1:n])
    v.bin <- names(validation.set[,test]) %in% names(v[1:n])
    c.table <- table(d.bin, v.bin)
    c.table['TRUE', 'TRUE']/n
  }) %>% unlist() %>%  enframe %>% 
    mutate(test=colnames(discovery.set))
  # AUROC-type?
  auc10 <- map(colnames(discovery.set), .f = function(test){
    d <- sort(discovery.set[,test])
    v <- validation.set[names(d), test]
    d[seq_along(d)] <- 0
    d[1:n] <- 1
    auc <- as.numeric(pROC::auc(predictor=v, response=d, 
                                levels=c(0,1), direction='>'))
  }) %>% unlist() %>%  enframe %>% 
    mutate(test=colnames(discovery.set))
  
  
  a <- map(colnames(discovery.set), .f = function(test){
    cor(-log10(discovery.set[,test]+1e-10), 
        -log10(validation.set[,test]+1e-10), 
        use='pairwise.complete.obs')
  }) %>% unlist() %>% 
    enframe %>% 
    mutate(test=colnames(discovery.set))
  res[[x]] <- full_join(c10, auc10 %>% rename(auc=value), by=c('name', 'test'))
}

res %>% 
  bind_rows() %>% 
  group_by(test) %>% 
  summarise(m=median(value), s=sd(value)) %>% 
  mutate(test=factor(test, levels=names(test.colours))) %>% 
  ggplot(aes(x=test, y=m, fill=test)) + 
    geom_bar(stat='identity') + 
    geom_point(data=res %>% bind_rows() %>%
                 mutate(test=factor(test, levels=names(test.colours))),
               aes(col=test, y=value)) +
    scale_fill_manual(values=alpha(test.colours, 0.5)) +
    # scale_fill_manual(values=test.colours) +
    scale_colour_manual(values=test.colours) +
    theme_publication() + 
    xlab('') + ylab('Median')
  
res %>% 
  bind_rows() %>% 
  group_by(test) %>% 
  summarise(m=median(auc), s=sd(auc)) %>% 
  mutate(test=factor(test, levels=names(test.colours))) %>% 
  ggplot(aes(x=test, y=m, fill=test)) + 
  geom_bar(stat='identity') + 
  geom_point(data=res %>% bind_rows() %>%
               mutate(test=factor(test, levels=names(test.colours))),
             aes(col=test, y=auc)) +
  scale_fill_manual(values=alpha(test.colours, 0.5)) +
  # scale_fill_manual(values=test.colours) +
  scale_colour_manual(values=test.colours) +
  theme_publication() + 
  xlab('') + ylab('AUC') + 
  coord_cartesian(ylim=c(0.5, 1))




# number of samples per study
n.all <- meta.crc %>% 
  group_by(Study) %>% 
  tally() %>% 
  bind_rows(meta.cd %>% 
    group_by(Study) %>% 
    tally())
# test confounder-aware
p.crc <- .f_test_naive_blocked(feat.crc, meta.crc, 'Group', 'Study')
true.crc <- p.crc$adj.blocked < 0.05
names(true.crc) <- p.crc$name
# p.cd <- .f_test_naive_blocked(feat.cd, meta.cd, 'Group', 'Study')
# true.cd <- p.cd$adj.blocked < 0.05
# names(true.cd) <- p.cd$name

# count number of significant
n.crc <- imap(all.results.crc, ~{
  p.fdr <- apply(.x[rownames(feat.crc),], 2, p.adjust, method='fdr')
  p.fdr[is.na(p.fdr)] <- 1
  vapply(colnames(p.fdr), function(x){
    pred <- p.fdr[,x] < 0.05
    if (sum(pred)==0){
      fpr <- 0
      tpr <- 0
    } else {
      tpr <- sum(names(which(pred)) %in% names(which(true.crc)))/sum(true.crc)
      fpr <- sum(!names(which(pred)) %in% names(which(true.crc)))/sum(!true.crc)
    }
    auroc <- as.numeric(pROC::auc(predictor=p.fdr[,x], response=true.crc,
                                  levels=c('FALSE', 'TRUE'), direction='<'))
    c('fpr'=fpr, 'tpr'=tpr, 'auroc'=auroc, 'n.discoveries'=sum(pred))
  }, FUN.VALUE = double(4)) %>% 
    t() %>% as_tibble(rownames = 'test') %>% 
    mutate(n=n.all %>% filter(Study%in% str_split(.y, pattern = '-')[[1]]) %>% 
             pull(n) %>% sum) %>% 
    mutate(studies=.y)}) %>% 
  bind_rows()
# n.cd <- imap(all.results.cd, ~{
#   p.fdr <- apply(.x, 2, p.adjust, method='fdr')
#   p.fdr[is.na(p.fdr)] <- 1
#   vapply(colnames(p.fdr), function(x){
#     pred <- p.fdr[,x] < 0.05
#     if (sum(pred)==0){
#       fpr <- 0
#       tpr <- 0
#     } else {
#       tpr <- sum(names(which(pred)) %in% names(which(true.cd)))/sum(true.cd)
#       fpr <- sum(!names(which(pred)) %in% names(which(true.cd)))/sum(!true.cd)
#     }
#     auroc <- as.numeric(pROC::auc(predictor=p.fdr[,x], response=true.cd,
#                                   levels=c('FALSE', 'TRUE'), direction='<'))
#     c('fpr'=fpr, 'tpr'=tpr, 'auroc'=auroc, 'n.discoveries'=sum(pred))
#   }, FUN.VALUE = double(4)) %>% 
#     t() %>% as_tibble(rownames = 'test') %>% 
#     mutate(n=n.all %>% filter(Study%in% str_split(.y, pattern = '-')[[1]]) %>% 
#              pull(n) %>% sum) %>% 
#     mutate(studies=.y)}) %>% 
#   bind_rows()

# plot N

g.n.crc <- n.crc %>% 
  filter(!str_detect(studies, '-')) %>% 
  filter(test!='ZIBSeq-sqrt') %>% 
  mutate(test=factor(test, levels = rev(names(test.colours)))) %>% 
  ggplot(aes(y=test, x=studies, fill=n.discoveries)) + 
  geom_tile() + 
  scale_fill_gradientn(colors=viridis::viridis(n=10)) + 
  ggembl::theme_publication() + 
  geom_text(aes(label=n.discoveries), colour='white') +
  xlab('') + ylab('')
ggsave(g.n.crc, filename = './figures/misc/real_data_results.pdf',
       width = 4, height = 4, useDingbats=FALSE)

n.discoveries <- n.crc %>% 
  filter(test!='ZIBSeq-sqrt') %>% 
  filter(n==max(n)) %>% 
  mutate(test=factor(test, levels=names(test.colours))) %>% 
  ggplot(aes(x=test, y=n.discoveries)) + 
    geom_bar(stat='identity', aes(fill=test)) + 
    scale_fill_manual(values=test.colours) +
    theme_publication()
ggsave(n.discoveries, filename = '~/Desktop/crc_all_test_n.pdf',
       width = 5, height = 4, useDingbats=FALSE)
n.all <- n.crc %>% 
  filter(test!='ZIBSeq-sqrt') %>% 
  # mutate(jp=str_detect(studies, 'Yachida')) %>% 
  mutate(test=factor(test, levels=names(test.colours))) %>% 
  ggplot(aes(x=n, y=n.discoveries, colour=test)) + 
    geom_point() + 
    geom_smooth(se=FALSE, method='glm')+
    scale_colour_manual(values=test.colours) +
    theme_publication()
ggsave(n.all, filename = '~/Desktop/crc_all_test_n_all.pdf',
       width = 6, height = 4, useDingbats=FALSE)
# plot TPR/FPR (against meta-analysis-type)
n.rates <- n.crc %>% 
  filter(test!='ZIBSeq-sqrt') %>% 
  # mutate(jp=str_detect(studies, 'Yachida')) %>% 
  filter(n==max(n)) %>% 
  mutate(test=factor(test, levels=names(test.colours))) %>% 
  ggplot(aes(y=tpr, x=fpr, col=test)) +
    geom_point(size=5) +
    scale_colour_manual(values=test.colours) + 
    theme_publication() + 
    xlab("False positive Rate") + 
    ylab("True positive Rate")
ggsave(n.rates, filename = '~/Desktop/crc_all_test_rates.pdf',
       width = 5, height = 4, useDingbats=FALSE)
auroc <- n.crc %>% 
  filter(test!='ZIBSeq-sqrt') %>% 
  mutate(test=factor(test, levels=names(test.colours))) %>% 
  mutate(auroc=1-auroc) %>% 
  mutate(jp=str_detect(studies, 'Yachida')) %>% 
  mutate(test=factor(test, levels=names(test.colours))) %>% 
  ggplot(aes(x=n, y=auroc, col=test)) +
    geom_point() +
    geom_smooth(se=FALSE) +
    scale_colour_manual(values=test.colours) + 
    theme_publication() + 
    facet_grid(~jp, scale='free_x')
ggsave(auroc, filename = '~/Desktop/crc_all_test_auroc.pdf',
       width = 6, height = 4, useDingbats=FALSE)

x <- all.results.crc$`Zeller_2014-Feng_2015-Yu_2017`
x[is.na(x)] <- 1
x.fdr <- apply(x, 2, p.adjust, method='fdr')
colSums(x.fdr < 0.05)

# TODO
# how to compare results?
  # number of discoveries per tool
  # similarities to confounder-aware tools?
  # overlap?

# ##############################################################################
# CRC stuff

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

# comparison of study/group effects in pairwise comparisons
pdf(here('figures', 'figure_batch_effects', 'real_comparisons_CRC.pdf'),
    width = 7, height = 3, useDingbats = FALSE)
studies <- unique(meta.crc$Study)
df.plot.phi <- tibble(Study1=character(0), Study2=character(0), phi=double(0))
for (i in seq_along(studies)){
  s1 <- studies[i]
  if (i==length(studies)) next()
  for (i2 in (i+1):length(studies)){
    s2 <- studies[i2]
    message(s1, '-', s2)
    meta.red <- meta.crc %>%  filter(Study%in%c(s1, s2))
    df.plot <- .f_contrast_gfc(feat.crc, meta.red,'Group', 'Study')
    df.test.w <- .f_test_naive_blocked(feat.crc, meta.red, 'Group', 
                                       'Study', 'wilcoxon')
    df.test.l <- .f_test_naive_blocked(feat.crc, meta.red, 'Group', 
                                       'Study', 'lme')
    df.test.limma <- .f_test_naive_blocked(feat.crc, meta.red, 
                                           'Group', 'Study', 'limma')
    g1 <- df.plot %>% 
      ggplot(aes(x=abs(label), y=abs(confounder))) + 
      geom_abline(slope = 1, intercept = 0) + 
      geom_point() + 
      xlim(0,2) + ylim(0,2) + 
      theme_publication() + 
      ggtitle(paste0(s1, '-', s2))
    df.test <- df.test.w %>% 
      mutate(adj.naive=p.adjust(naive, method='fdr')) %>% 
      mutate(adj.blocked=p.adjust(blocked, method='fdr')) %>% 
      mutate(found.naive=adj.naive < 0.05, found.blocked=adj.blocked < 0.05) %>% 
      mutate(test='Wilcoxon') %>% 
      bind_rows(df.test.l %>% 
      mutate(adj.naive=p.adjust(naive, method='fdr')) %>% 
      mutate(adj.blocked=p.adjust(blocked, method='fdr')) %>% 
      mutate(found.naive=adj.naive < 0.05, 
             found.blocked=adj.blocked < 0.05) %>% mutate(test='LM')) %>% 
      bind_rows(df.test.limma %>% 
      mutate(adj.naive=p.adjust(naive, method='fdr')) %>% 
      mutate(adj.blocked=p.adjust(blocked, method='fdr')) %>% 
      mutate(found.naive=adj.naive < 0.05, 
             found.blocked=adj.blocked < 0.05) %>% mutate(test='limma'))
    
    sp <- df.test %>% 
      group_by(test) %>% 
      summarise(m=cor(blocked, naive, method='spearman', 
                      use='pairwise.complete.obs'))
    g2 <- df.test %>% 
      ggplot(aes(x=-log10(naive), y=-log10(blocked), col=test)) +
      geom_abline(slope = 1, intercept = 0) + 
      geom_point() + 
      theme_publication() +
      annotate(geom='text', x=5, y=c(8, 10, 12),
               label=paste0('rho_', sp$test, '=', sprintf(fmt='%.2f', sp$m)))
    g <- plot_grid(g1, g2, nrow = 1)
    print(g)
    
    df.plot.phi <- df.plot.phi %>%
      add_row(Study1=s1, Study2=s2,
              phi=-phi(table(meta.red$Group, meta.red$Study)))
  }
}
dev.off()

  
# ##############################################################################
# CD-stuff

# convert to genus level abundances
feat.cd.genus <- read.table('data/motus_cd_genus.tsv', sep='\t', header=TRUE,
                            check.names = FALSE)
feat.cd.genus <- prop.table(as.matrix(feat.cd.genus), 2)

pdf(here('figures', 'figure_batch_effects', 'real_comparisons_CD.pdf'),
    width = 7, height = 3, useDingbats = FALSE)
studies <- unique(meta.cd$Study)
for (i in seq_along(studies)){
  s1 <- studies[i]
  if (i==length(studies)) next()
  for (i2 in (i+1):length(studies)){
    for (tbl in c('motus', 'genus')){
      if (tbl == 'genus'){
        f.cd <- feat.cd.genus
      } else {
        f.cd <- feat.cd
      }
    s2 <- studies[i2]
    meta.red <- meta.cd %>% 
      filter(Study%in%c(s1, s2))
    df.plot <- .f_contrast_gfc(f.cd, meta.red, 'Group', 'Study')
    df.test.w <- .f_test_naive_blocked(f.cd, meta.red, 'Group', 
                                       'Study', 'wilcoxon')
    df.test.l <- .f_test_naive_blocked(f.cd, meta.red, 'Group', 
                                       'Study', 'lme')
    df.test.limma <- .f_test_naive_blocked(f.cd, meta.red, 'Group', 
                                           'Study', 'limma')
    g1 <- df.plot %>% 
      ggplot(aes(x=abs(label), y=abs(confounder))) + 
      geom_abline(slope = 1, intercept = 0) + 
      geom_point() + 
      xlim(0,2) + ylim(0,2) + 
      theme_publication() + 
      ggtitle(paste0(s1, '-', s2))
    df.test <- df.test.w %>% 
      mutate(adj.naive=p.adjust(naive, method='fdr')) %>% 
      mutate(adj.blocked=p.adjust(blocked, method='fdr')) %>% 
      mutate(found.naive=adj.naive < 0.05, found.blocked=adj.blocked < 0.05) %>% 
      mutate(test='Wilcoxon') %>% 
      bind_rows(df.test.l %>% 
                  mutate(adj.naive=p.adjust(naive, method='fdr')) %>% 
                  mutate(adj.blocked=p.adjust(blocked, method='fdr')) %>% 
                  mutate(found.naive=adj.naive < 0.05, 
                         found.blocked=adj.blocked < 0.05) %>% 
                  mutate(test='LM')) %>% 
      bind_rows(df.test.limma %>% 
                  mutate(adj.naive=p.adjust(naive, method='fdr')) %>% 
                  mutate(adj.blocked=p.adjust(blocked, method='fdr')) %>% 
                  mutate(found.naive=adj.naive < 0.05, 
                         found.blocked=adj.blocked < 0.05) %>% 
                  mutate(test='limma'))
    
    df.test %>% 
      select(name, found.naive, found.blocked, test) %>% 
      full_join(df.plot, by='name') %>% 
      mutate(type=case_when(found.naive & found.blocked ~ 'unaffected',
                            found.naive & !found.blocked ~ 'affected',
                            !found.naive & found.blocked ~ 'hidden',
                            TRUE~'background')) %>% 
      ggplot(aes(x=abs(label), y=abs(confounder), col=type)) + 
      geom_abline(slope = 1, intercept = 0) + 
      geom_point()  +
        facet_grid(~test)
    
    g <- df.test %>% 
      mutate(type=case_when(found.naive & found.blocked ~ 'unaffected',
                            found.naive & !found.blocked ~ 'affected',
                            !found.naive & found.blocked ~ 'hidden',
                            TRUE~'background')) %>% 
      ggplot(aes(x=-log10(adj.naive), -log10(adj.blocked), col=type)) + 
      geom_abline(slope = 1, intercept = 0) + 
        geom_hline(yintercept = -log10(0.05)) +
        geom_vline(xintercept = -log10(0.05)) +
        geom_point() + 
        facet_grid(~test) + 
        theme_publication(panel.grid = 'major')
    
    
    
    sp <- df.test %>% 
      group_by(test) %>% 
      summarise(m=cor(blocked, naive, method='spearman', 
                      use='pairwise.complete.obs'))
    g2 <- df.test %>% 
      ggplot(aes(x=-log10(naive), y=-log10(blocked), col=test)) +
      geom_abline(slope = 1, intercept = 0) + 
      geom_point() + 
      theme_publication() +
      annotate(geom='text', x=5, y=c(15, 16, 17),
               label=paste0('rho_', sp$test, '=', sprintf(fmt='%.2f', sp$m)))
    g <- plot_grid(g1, g2, nrow = 1)
    print(g)
    }
    df.plot.phi <- df.plot.phi %>%
      add_row(Study1=s1, Study2=s2,
              phi=-phi(table(meta.red$Group, meta.red$Study)))
  }
}
dev.off()
