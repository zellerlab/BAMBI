# ##############################################################################
#
## Assess the effect sizes of true confounders
##
##  Alcohol in AGP
##  Metformin in T2D
##  Metacardis?
##  Abx intake in Milieu Interieur?
#
# ##############################################################################

library("here")
library("tidyverse")
library("vegan")


.f_conf_tbls <- function(meta, feat, sampleid, label){
  
  # tests
  stopifnot(sampleid %in% colnames(meta))
  if (!is.na(label)) stopifnot(label %in% colnames(meta))
  stopifnot(all(meta[[sampleid]] %in% colnames(feat)))
  feat <- feat[,meta[[sampleid]]]
  # assess confounder bias --> wilcox.test/fisher.test (like SIAMCAT)
  # assess global effect of confounders --> adonis-test
  
  temp <- prop.table(as.matrix(feat), 2)
  temp <- temp[rowMeans(temp!=0) > 0.05,]
  # temp.dist <- vegan::vegdist(t(temp),method='bray')
  
  df.meta <- tibble(conf=character(0), p.adonis=double(0), f.adonis=double(0))
  message('+++ calculating global confounder influence')
  for (x in colnames(meta)){
    message(x)
    if (x == sampleid){
      next()
    }
    conf <- meta[[x]]
    names(conf) <- meta[[sampleid]]
    if (any(is.na(conf))){
      conf <- conf[!is.na(conf)]
    }
    res.adonis <- vegan::adonis(
      t(temp[,names(conf)])~as.factor(conf))
    df.meta <- df.meta %>% 
      add_row(conf=x, 
              p.adonis=res.adonis$aov.tab$`Pr(>F)`[1],
              f.adonis=res.adonis$aov.tab$F.Model[1])
  }
  if (!is.na(label)){
    message('+++ calculating associations between confounders and label')
    p.val <- tibble(conf=character(0), p.value=double(0))
    for (x in colnames(meta)){
      if (x %in% c(sampleid, label)){
        next()
      }
      if (sum(is.na(meta[[x]])) > 100){
        p.val <- p.val %>% 
          add_row(conf=x, p.value=NA_real_)
        next()
      }
      if (typeof(meta[[x]]) == 'character'){
        if (length(unique(na.omit(meta[[x]]))) > 1){
          if (length(unique(meta[[x]])) > 10){
            t <- list(p.value=NA)
          } else {
            t <- fisher.test(meta[[label]], meta[[x]])
          }
        } else {
          t <- list(p.value=NA)
        }
      } else if (typeof(meta[[x]]) == 'double'){
        t <- wilcox.test(meta[[x]]~meta[[label]])
      }
      p.val <- p.val %>% 
        add_row(conf=x, p.value=t$p.value)
    }
    df.meta <- df.meta %>% full_join(p.val, by='conf')
  }
  
  # effect size for multiple group confounders --> variance associated with conf
  df.var <- list()
  message('+++ calculating associations between confounders and features')
  for (x in colnames(meta)){
    message(x)
    if (x==sampleid){
      next()
    }
    conf <- meta[[x]]
    names(conf) <- meta[[sampleid]]
    if (any(is.na(conf))){
      conf <- conf[!is.na(conf)]
    }
    if (x=="BMI"){
      conf.d <- conf
      conf <- case_when(conf.d < 18.5~'U',
                        conf.d >= 18.5 & conf.d <= 24.9 ~ "H",
                        conf.d > 24.9 & conf.d <= 29.9 ~ "Ov",
                        conf.d >29.9 ~ "Ob")
      names(conf) <- names(conf.d)
    } else if (is.double(conf)){
      quart <- quantile(conf, probs = seq(0, 1, 0.25), na.rm = TRUE)
      conf.q <- cut(conf, unique(quart), include.lowest = TRUE)
      names(conf.q) <- names(conf)
      conf <- conf.q
    }
    var.all <- vapply(rownames(temp), FUN=function(f){
      x.rank <- temp[f,names(conf)]
      x.rank <- rank(x.rank)/length(x.rank)
      ss.tot <- sum((x.rank - mean(x.rank))^2)/length(x.rank)
      ss.o.i <- sum(vapply(unique(conf), function(s){
        sum((x.rank[conf==s] - mean(x.rank[conf==s]))^2)
      }, FUN.VALUE = double(1)))/length(x.rank)
      return(1-ss.o.i/ss.tot)
    }, FUN.VALUE = double(1))
    df.var[[x]] <- var.all
  }
  
  # plot one against the other
  df.plot <- bind_cols(df.var, .id=rownames(temp))
  
  return(list('conf'=df.meta, 'var'=df.plot))
}

.f_rsq <- function(actual, pred){
  rss <- sum((pred - actual) ^ 2)  ## residual sum of squares
  tss <- sum((actual - mean(actual)) ^ 2)  ## total sum of squares
  rsq <- 1 - rss/tss
  return(rsq)
}


# ##############################################################################
# AGP data
meta.agp <- read_tsv(here('data', 'meta_AGP_2020.tsv'))
feat.agp <- read.table(here('data', 'otus_AGP_2020.tsv'), 
                       sep='\t', stringsAsFactors = FALSE, check.names = FALSE,
                       quote = '', comment.char = '')
feat.agp <- as.matrix(feat.agp)
meta.agp.red <- meta.agp %>% 
  mutate(country=case_when(country_USA==1~'USA',
                           country_Canada==1~'Canada',
                           country_United.Kingdom==1~'UK',
                           TRUE~NA_character_)) %>% 
  select(sample_name, bmi, age_years, sex, alcohol_frequency, country,
         bowel_movement_quality, country,
         milk_cheese_frequency, vegetable_frequency, antibiotic_history)


# robustness of confounder-associations
# log rel.ab features
# LM for associations with meta-variables
# check residuals on left-out data?
feat <- prop.table(feat.agp, 2)
feat <- log10(feat + 1e-05)

combined.res <- list()

for (i in seq_len(rep)){
  # subsample dataset
  meta.test <- meta.agp.red %>% 
    slice_sample(n=nrow(meta.agp.red)*0.4)
  meta.train <- meta.agp.red %>% 
    filter(!sample_name %in% meta.test$sample_name)
  
  df.train <- t(feat) %>% as_tibble(rownames = 'sample_name') %>% 
    right_join(meta.train, by='sample_name') %>% 
    as.data.frame()
  colnames(df.train) <- make.names(colnames(df.train))
  
  df.test <- t(feat) %>% as_tibble(rownames = 'sample_name') %>% 
    right_join(meta.test, by='sample_name') %>% 
    as.data.frame()
  colnames(df.test) <- make.names(colnames(df.test))
  
  # loop through meta-variables
  for (var in setdiff(colnames(meta.agp.red), 'sample_name')){
    message(var)
    
    # loop through features
    if (any(is.na(df.train[[var]]))){
      df.train.red <- df.train[!is.na(df.train[[var]]),]
    } else {
      df.train.red <- df.train
    }
    
    if (any(is.na(df.test[[var]]))){
      df.test.red <- df.test[!is.na(df.test[[var]]),]
    } else {
      df.test.red <- df.test
    }
    
    # train linear model
    models <- lapply(rownames(feat), FUN=function(x){
      x <- make.names(x)
      f <- paste0(x, '~', var)
      fit <- lm(formula=f, data=df.train.red)})
    names(models) <- make.names(rownames(feat))
    rsq.train <- map_dbl(models, .f=function(x){summary(x)$r.squared})
    
    
    # apply to test set
    rsq.test <- map_dbl(names(models), .f = function(x){
      preds <- predict(models[[x]], newdata=df.test.red)
      .f_rsq(actual = df.test[[x]], pred = preds)
    })
    
    # record 
    combined.res[[length(combined.res) + 1]] <- 
      tibble(trainig=rsq.train, test=rsq.test, variable=var, resample=i,
             feat.id=names(rsq.train))
  }
}
df.plot <- combined.res %>% bind_rows()

















# take the T2D cohorts
# unmatched
meta.agp.t2d.um <- read_csv('./figures/diabetes_unmatched.csv') %>% 
  mutate_all(as.character) %>% 
  mutate(bmi=as.numeric(bmi), age_years=as.numeric(age_years)) %>% 
  select(-longitude, -latitude, -pairID, -pairDist, -worstPairDist)
# matchted
meta.agp.t2d.m <- read_csv('./figures/diabetes_matched.csv') %>% 
  mutate_all(as.character) %>% 
  mutate(bmi=as.numeric(bmi), age_years=as.numeric(age_years)) %>% 
  select(-longitude, -latitude, -pairID, -pairDist, -worstPairDist)

res.t2d.m <- .f_conf_tbls(meta.agp.t2d.m, feat.agp, 'sample_name', 'target')
res.t2d.um <- .f_conf_tbls(meta.agp.t2d.um, feat.agp, 'sample_name', 'target')

res.t2d.um$var %>% 
  pivot_longer(-c(.id, target), names_to = 'conf', values_to = 'variance') %>% 
  ggplot(aes(x=target, y=variance)) + 
    geom_abline(slope = 1, intercept = 0) +
    geom_point() + 
    facet_wrap(~conf)
res.t2d.m$var %>% 
  pivot_longer(-c(.id, target), names_to = 'conf', values_to = 'variance') %>% 
  ggplot(aes(x=target, y=variance)) + 
  geom_point() + 
  facet_wrap(~conf)


x <- res.t2d.m$var %>% arrange(desc(alcohol_frequency)) %>% slice(3) %>% pull(.id) 
prop.table(feat.agp, 2)[x,meta.agp.t2d.m$sample_name] %>% 
  enframe(name='sample_name', value = 'feat') %>% 
  left_join(meta.agp.t2d.m) %>% 
  ggplot(aes(x=alcohol_frequency, y=log10(feat+1e-04))) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.1)


g1 <- res.t2d.um$var %>% 
  pivot_longer(-.id, names_to='conf', values_to='var.unmatched') %>% 
  full_join(res.t2d.m$var %>% 
              pivot_longer(-.id, names_to='conf', values_to='var.matched')) %>% 
  pivot_longer(-c(.id, conf)) %>% 
  ggplot(aes(x=conf, y=value, fill=name)) + 
    geom_boxplot() + 
    coord_flip() + 
    scale_fill_manual(values=c('#009F4D', '#307FE2'), guide=FALSE)

g2 <- res.t2d.um$conf %>% 
  transmute(conf, f.um=f.adonis) %>% 
  full_join(res.t2d.m$conf %>% 
              transmute(conf, f.m=f.adonis)) %>% 
  pivot_longer(-conf) %>% 
  ggplot(aes(x=conf, y=value, fill=name)) + 
    geom_bar(stat='identity', position = position_dodge()) + 
    coord_flip() + 
    scale_fill_manual(values=c('#009F4D', '#307FE2')) + 
    theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())
cowplot::plot_grid(g1, g2)


res.t2d.um$var %>% 
  transmute(bmi.um=alcohol_frequency, .id) %>% 
  full_join(
    res.t2d.m$var %>% 
      transmute(bmi.m=alcohol_frequency, .id)) %>% 
  ggplot(aes(x=bmi.um, y=bmi.m)) + 
    geom_point()



# ##############################################################################
# check overall in the AGP data
meta.agp.red <- meta.agp %>% 
  select(sample_name, sex, age_years, bmi, dog, nail_biter, 
         antibiotic_history, 
         bowel_movement_frequency, bowel_movement_quality, 
         alcohol_consumption, alcohol_frequency, 
         liver_disease, diabetes, ibd, cardiovascular_disease, cdiff, 
         asd, sibo, ibs, cancer, acid_reflux,
         diet_type_vegetarian, diet_type_vegam, diet_type_pescatarian, 
         diet_type_omnivore, diet_type_omnivore_nored, 
         country_Australia, country_Canada, country_USA, country_United.Kingdom,
         collection_season_Fall, collection_season_Spring, 
         collection_season_Summer, collection_season_Winter) %>% 
  mutate_all(as.character) %>% 
  mutate(bmi=as.numeric(bmi), age_years=as.numeric(age_years)) %>% 
  mutate(collection_season=case_when(collection_season_Spring=='1'~'Spring',
                                     collection_season_Summer=='1'~'Summer',
                                     collection_season_Fall=='1'~'Fall',
                                     collection_season_Winter=='1'~'Winter')) %>% 
  mutate(country=case_when(country_Canada=='1'~'Canada',
                           country_USA=='1'~'USA',
                           country_Australia=='1'~'Australia',
                           country_United.Kingdom=='1'~'UK')) %>% 
  mutate(diet=case_when(diet_type_vegetarian=='1'~'vegetarian',
                        diet_type_vegam=='1'~'vegan',
                        diet_type_pescatarian=='1'~'pescatarian',
                        diet_type_omnivore=='1'~'omnivore',
                        diet_type_omnivore_nored=='1'~'omnivore_nored')) %>% 
  select(-c(diet_type_vegetarian, diet_type_vegam, diet_type_pescatarian, 
            diet_type_omnivore, diet_type_omnivore_nored, 
            country_Australia, country_Canada, country_USA, 
            country_United.Kingdom,
            collection_season_Fall, collection_season_Spring, 
            collection_season_Summer, collection_season_Winter))
# res.agp.all <- .f_conf_tbls(meta.agp.red, feat.agp, 
#                             'sample_name', NA_character_)
# df.plot <- res.agp.all$var
# 
# df.plot %>% 
#   pivot_longer(-.id) %>% 
#   ggplot(aes(x=name, y=value)) + 
#     geom_boxplot() + 
#     coord_flip()
# 
# x <- df.plot %>% arrange(desc(alcohol_frequency)) %>% slice(1) %>% pull(.id)
# meta.agp.red %>% 
#   left_join(enframe(prop.table(feat.agp, 2)[x,], 
#                     name='sample_name', value='feat')) %>% 
#   filter(!is.na(alcohol_frequency)) %>% 
#   ggplot(aes(x=alcohol_frequency, y=log10(feat+1e-04))) + 
#     geom_boxplot(outlier.shape = NA) + 
#     geom_jitter(width = 0.1)




# ##############################################################################
# Metacardis data
meta.mc <- read_tsv(here('data', 'meta_t2d_metacardis.tsv'))
feat.mc <- read.table(here('data', 'motus_t2d_metacardis.tsv'), 
                      quote = '', comment.char = ',', row.names = 1,
                      stringsAsFactors = FALSE, check.names = FALSE, 
                      header = TRUE)
meta.mc.red <- meta.mc %>% 
  select(SampleID, LABEL, AGE, ANTIBIOTIC, BMI_C, STATINE_C, METFORMIN_C,
         SMOKSTATUS_C, PPI_C, nutr_alcohol)
res.mc <- .f_conf_tbls(meta.mc.red, feat.mc, 'SampleID', 'LABEL')
df.plot <- res.mc$var

df.plot %>% 
  pivot_longer(-c(.id, LABEL)) %>% 
  ggplot(aes(x=LABEL, y=value)) + 
    geom_point() +
    facet_wrap(~name)

# ##############################################################################
# Metformin data
meta.metf <- read_tsv(here('data', 'meta_metformin.tsv'))
feat.metf <- read.table(here('data', 'motus_metformin.tsv'),
                        sep='\t', stringsAsFactors = FALSE, check.names = FALSE,
                        quote = '', comment.char = '')
feat.metf <- as.matrix(feat.metf)

res.metf <- .f_conf_tbls(meta.metf, feat.metf, 'Sample_ID', 'GROUP')
df.plot <- res.metf$var

df.plot %>% 
  # select(-COUNTRY) %>% 
  pivot_longer(-c(.id, GROUP)) %>% 
  ggplot(aes(x=GROUP, y=value)) + 
    geom_point() +
    facet_wrap(~name)



# ##############################################################################
# mouse data from the mouse gene catalgoue?

meta.mice <- read_tsv('./data/meta_Xiao_2015.tsv', col_names = FALSE)
feat.mice <- read.table(here('data', 'Xiao_2015.motus'),
                        sep='\t', stringsAsFactors = FALSE, check.names = FALSE,
                        quote = '', comment.char = '', header = 1, 
                        row.names = 185)
feat.mice <- feat.mice[,colSums(feat.mice) > 1000]
feat.mice <- as.matrix(feat.mice)

# select relevant metadata
# age, diet, supplier, country, sex, genetics
meta.mice.red <- meta.mice %>% 
  select(X6, X10, X15, X16, X20, X21, X22, X23) %>% 
  rename(Age=X6, Diet=X10, Food_supplier=X15, Country=X16,Sex=X20,
         Lineage=X21, Mass=X22, Sample_ID=X23) %>% 
  mutate(Age=as.numeric(str_remove(Age, 'Age:'))) %>%
  mutate(Diet=str_remove(Diet, 'Diet:')) %>% 
  mutate(Food_supplier=str_remove(Food_supplier, 'Food supplier:')) %>% 
  mutate(Sex=str_remove(Sex, 'Sex:')) %>% 
  mutate(Lineage=str_remove(Lineage, "Subspecific genetic lineage:")) %>% 
  mutate(Country=str_remove(Country, "Geographic location \\(country and/or sea,region\\):")) %>% 
  mutate(Mass=as.numeric(str_remove(Mass, 'Total mass:'))) %>% 
  mutate(Sample_ID=str_remove(Sample_ID, 'Alternative accession-BioSample:')) %>% 
  filter(Sample_ID %in% colnames(feat.mice))

res.mice <- .f_conf_tbls(meta.mice.red, feat.mice, 'Sample_ID', 'Diet')
res.mice$var %>% 
  pivot_longer(-c(.id, Diet)) %>% 
  ggplot(aes(x=Diet, y=value)) + 
    geom_point() + 
    facet_wrap(~name)
  
id <- 'Rikenellaceae species incertae sedis [ext_mOTU_v26_18923]'
meta.mice.red %>% 
  mutate(feat=prop.table(feat.mice, 2)[id, Sample_ID]) %>% 
  mutate(feat=log10(feat+1e-05)) %>% 
  ggplot(aes(x=Country, y=feat, fill=Country)) + 
    geom_boxplot()
