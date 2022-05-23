library("here")
library("tidyverse")
library("ggembl")

# load results (morgan)
# lm, lme, wilcoxon, ANCOMBC -- jakob's with all norms
# _ext -- with rank, rarefy, rarefy-TSS norms
# _ext -- with pass and TSS.log for metadeconfoundR lmes
df.tmp <- map_dfr(list.files(here('test_results','conf_study'), 
                             pattern = 'conf_study'),
                  ~ read_tsv(here('test_results', 'conf_study', .))) %>%
  mutate(norm = str_replace_all(norm, 'rarefy-TSS', 'rarefy.TSS'))

# filter to match j's specifications (morgan)
df.red <- df.tmp %>% 
  filter(str_detect(problem, paste(c('0.5','0.7','0.9'),
                                   collapse = '|'))) %>%
  filter(str_detect(group, 'ab4_prev3')) %>% 
  filter(subset==200) %>% 
  select(-c(mem.used, time.running)) %>% 
  mutate(combi=paste0(test, '-', norm)) %>% 
  filter(combi %in% c("limma-pass", 'limma_conf-pass',
                      "wilcoxon-pass", "wilcoxon_conf-pass",
                      "wilcoxon-TSS.log", "wilcoxon_conf-TSS.log",
                      "wilcoxon-rarefy", "wilcoxon_conf-rarefy",
                      # wouldn't this lme behave more similarly to mdc?
                      "lm-rank", "lme_conf-rank",
                      "lm-TSS.log", "lme_conf-TSS.log",
                      #"lm-rarefy.TSS","lme_conf-rarefy.TSS",
                      "mdc-FE_conf-pass", "mdc-RE_conf-pass",
                      "mdc-FE_conf-TSS.log", "mdc-RE_conf-TSS.log",
                      "mdc-FE_conf-rarefy.TSS", "mdc-RE_conf-rarefy.TSS",
                      "ANCOMBC-pass", "ANCOMBC_conf-pass")) %>%
  mutate(SENS = TP/(TP+FN))
write_tsv(df.red, here('figures', 'figure_conf_study', 
                            'test_results_red.tsv')) 

method.ranking.all <- df.tmp %>%
  filter(str_detect(problem, paste(c('0.5','0.7','0.9'),
                                   collapse = '|'))) %>%
  separate(group, into=c('group', 'repetition'), sep='_rep') %>%
  filter(group!='ab1_prev1') %>% 
  mutate(subset=as.numeric(as.character(subset))) %>% 
  filter(subset > 30 & subset < 300) %>% 
  group_by(test, norm) %>% 
  summarise(n.10=sum(FDR > 0.1), m.auroc=mean(AUC), .groups = 'drop') %>% 
  arrange(n.10, desc(m.auroc)) %>% 
  group_by(test) %>% 
  slice(n=1) %>% 
  mutate(type=paste0(test, '-', norm))

df.comb <- df.tmp %>%
  mutate(PR = 1-FDR) %>%
  select(-job.id, -time.running, -mem.used) %>%
  mutate(combi = paste0(test, '-', norm)) %>% 
  filter(combi %in% c("limma-pass", 'limma_conf-pass',
                      "wilcoxon-pass", "wilcoxon_conf-pass",
                      "wilcoxon-TSS.log", "wilcoxon_conf-TSS.log",
                      "wilcoxon-rarefy", "wilcoxon_conf-rarefy",
                      # wouldn't this lme behave more similarly to mdc?
                      "lm-rank", "lme_conf-rank",
                      "lm-TSS.log", "lme_conf-TSS.log",
                      #"lm-rarefy.TSS","lme_conf-rarefy.TSS",
                      "mdc-FE_conf-pass", "mdc-RE_conf-pass",
                      "mdc-FE_conf-TSS.log", "mdc-RE_conf-TSS.log",
                      "mdc-FE_conf-rarefy.TSS", "mdc-RE_conf-rarefy.TSS",
                      "ANCOMBC-pass", "ANCOMBC_conf-pass")) %>%
  mutate(SENS = TP/(TP+FN)) %>%
  mutate(conf = case_when(str_detect(test, fixed('_conf')) ~ 'corrected', 
                          TRUE ~ 'naive')) %>% 
  mutate(conf = factor(conf, levels = c('naive', 'corrected'))) %>% 
  mutate(test = str_remove(test, '_conf')) %>% 
  mutate(combi = str_remove(combi, '_conf')) %>% 
  mutate(test = case_when(test %in% c('lme',
                                      'mdc-FE',
                                      'mdc-RE') ~ 'lm', 
                          TRUE ~ test)) %>% 
  group_by(problem, group, subset, test, combi, conf) %>%
  summarise(precision=mean(PR), sd.precision=sd(PR),
            recall=mean(R), sd.recall=sd(R),
            AUC=mean(auroc), sd.AUC=sd(auroc),
            FDR=mean(FDR), sd.fdr=sd(FDR), SENS=mean(SENS),
            .groups='drop') %>%
  separate(group, into=c('ab','prev'), sep = '_', remove=FALSE) 
write_tsv(df.comb,
          file = here('figures', 'figure_conf_study',
                      'test_results_ext.tsv'))

# df.comb <- read_tsv(here('figures', 'figure_conf_study',
#                          'test_results_ext.tsv'))
# df.red <- read_tsv(here('figures', 'figure_conf_study', 
#                         'test_results_red_ext.tsv'))

df.plot <- df.red %>%
  mutate(conf = case_when(str_detect(test, fixed('_conf')) ~ 'corrected', 
                          TRUE ~ 'naive')) %>% 
  mutate(conf = factor(conf, levels = c('naive', 'corrected'))) %>% 
  mutate(test = str_remove(test, '_conf')) %>% 
  mutate(combi = str_remove(combi, '_conf')) %>% 
  mutate(test = case_when(test %in% c('lme',
                                      'mdc-FE',
                                      'mdc-RE') ~ 'lm', 
                          TRUE ~ test)) %>% 
  select(problem, test, SENS, FDR, conf, combi) %>% 
  # mutate(auroc = (auroc-0.5)/0.5) %>% 
  pivot_longer(cols = c(SENS, FDR)) #%>%
  #filter(!str_detect(combi, 'mdc-RE'))

plot.fdr <- df.plot %>%
  filter(name == 'FDR') %>%
  ggplot(aes(x=combi, y=value, alpha=conf, fill='red')) + 
  geom_boxplot() + 
  theme_bw() +
  labs(x = '', y = '') +
  guides(alpha = 'none', fill = 'none',title='FDR') +
  facet_grid(rows = vars(problem), cols = vars(test), scales = 'free') +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 9),
        axis.line.x = element_blank())

plot.sens <- df.plot %>%
  filter(name == 'SENS') %>%
  ggplot(aes(x=combi, y=value, alpha=conf, fill = 'blue')) + 
  geom_boxplot() + 
  scale_fill_manual(values = 'blue') +
  theme_bw() +
  labs(x = '', y = '',title='Sensitivity') +
  guides(alpha = 'none',fill = 'none') +
  facet_grid(rows = vars(problem), cols = vars(test), scales = 'free') +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 9),
        axis.line.x = element_blank())
  
plot.both <- df.plot %>%
  #filter(test %in% c('lm','wilcoxon')) %>%
  ggplot(aes(x=combi, y=value, alpha=conf, fill=name)) + 
  geom_hline(yintercept = 0.5, color = 'grey') +
  geom_boxplot() + 
  labs(x = '', y = '') +
  scale_alpha_manual(values=c(0.2, 1)) + 
  scale_fill_manual(values=c('#E40046','#307FE2')) +
  theme_bw() +
  labs(fill = 'Measurement') + guides(alpha = 'none') +
  facet_grid(rows = vars(problem), cols = vars(test), scales = 'free')+
  theme(axis.text.x = element_text(angle = -45, hjust = 0, size = 9),
          axis.line.x = element_blank())
