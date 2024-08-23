# ##############################################################################
#
## Supplementary Figure: Reality assessment and vanilla results for 
##    compositionally-confounded simulations
#
# ##############################################################################

library("tidyverse")
library("here")
library("ggembl")
library("ggrepel")

# ##############################################################################
# effect size plot

effect.size.all <- read_tsv(here('reality_checks', 'sim_Zeevi_WGS_all', 
                                 'variance.tsv'))
effect.size.comp <- read_tsv(here('reality_checks', 
                                  'sim_Zeevi_WGS_compositional', 
                                  'variance.tsv'))

# part A: shift in background feature effect size
g <- effect.size.all %>% 
  bind_rows(effect.size.comp) %>% 
  filter(prev==2, rep==89) %>% 
  ggplot(aes(x=as.factor(ab), y=fc.sim, fill=selection)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1),
                size=0.7, aes(col=selection)) + 
    facet_grid(~simulation) + 
    theme_publication(panel.grid = 'major_y') + 
    xlab('Abundance scaling') + 
    ylab('Generalized fold change') + 
    scale_fill_manual(values=alpha(c('#707372', '#A8C700'), alpha=0.4), 
                      name='Implanted') +
    scale_colour_manual(values=c('#707372', '#A8C700'), name='Implanted')
ggsave(g, filename = here('figures', 'figure_comp', 'effect_size.pdf'),
       width = 6, height = 4, useDingbats=FALSE)

# ##############################################################################
# results from vanilla testing

# take a specific group/sample size
fn.final <- here('figures', 'figure_comp', 'test_results.tsv')
if (!file.exists(fn.final)){
  .f_read_single_test <- function(fn, fn.final){
    message(fn)
    res <- vroom::vroom(fn) %>% 
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
  
  fn.test.results <- list.files(here('test_results'), 
                                pattern='compositional.*.tsv',
                                full.names = TRUE, recursive = TRUE)
  df.res.all <- map(fn.test.results, .f=.f_read_single_test, fn.final=fn.final)
} else {
  df.res.all <- read_tsv(fn.final)
}

test.colours <- yaml::yaml.load_file(here('files', 'test_colours.yml'))
params <- yaml::yaml.load(here('create_simulations', 'parameters.yaml'))

g <- df.res.all %>% 
  bind_rows(df.res.all %>% 
              filter(test=='ANCOM') %>% 
              mutate(adjust='BH')) %>% 
  filter(prev=='prev2', subset==200) %>% 
  mutate(g=paste0(test, '-', norm)) %>% 
  filter(g %in% c('ALDEx2-pass', 'ANCOM-rarefy', 'ANCOMBC-pass', 'corncob-pass', 
                  'DESeq2-rarefy', 'distinct-pass', 'edgeR-pass', 'fastANCOM-pass', 
                  'KS-clr', 'LDM-pass', 'limma-TSS.log', 'LinDA-pass', 'lm-TSS.log', 
                  'metagenomeSeq-pass', 'metagenomeSeq2-pass', 'wilcoxon-TSS', 
                  'ZIBSeq-sqrt', 'ZicoSeq-pass', 'ZINQ-pass', 
                  't-test-TSS.log')) %>% 
  dplyr::select(ab, test, precision, recall, AUC, fpr, adjust) %>% 
  mutate(test=factor(test, levels = names(test.colours))) %>% 
  mutate(fpr=1-fpr) %>% 
  pivot_longer(c(precision, recall, AUC, fpr)) %>% 
  ggplot(aes(x=ab, y=value)) + 
    geom_line(aes(group=test, colour=test)) +
    scale_colour_manual(values=unlist(test.colours)) +
    facet_grid(adjust~name, scales = 'free_y') + 
    theme_presentation(panel.grid = 'major_y')
    
ggsave(g, filename = here('figures', 'figure_comp', 'all_measures.pdf'),
       width = 5, height = 4, useDingbats=FALSE)
