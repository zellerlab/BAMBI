# ##############################################################################
#
## Run tests for real CRC and CD data
#
# ##############################################################################

library("tidyverse")
library("batchtools")
library("here")

# load SIMBA
args <- commandArgs(trailingOnly = TRUE)
disease <- args[1]
stopifnot(disease %in% c('CRC', 'CD'))
devtools::load_all(simba.loc)

tests <- c('ANCOM', 'ANCOMBC', 'corncob', 'ZIBSeq', 'ZIBSeq-sqrt', 'DESeq2',
           'edgeR', 'metagenomeSeq', 'metagenomeSeq2', 'lm', 'wilcoxon',
           'limma', 'distinct', 'ZINQ', 'KS', 'ALDEx2')
if (!dir.exists(paste0(temp.loc, 'test_results_registries/'))){
  dir.create(paste0(temp.loc, 'test_results_registries/'))
}
job.registry <- paste0(temp.loc, 'test_results_registries/real_data_', disease)
if (!dir.exists(job.registry)) {
  reg <- makeRegistry(file.dir = job.registry, 
                                work.dir = here(),
                                conf.file = here('cluster_config',
                                                 'batchtools_test_conf.R'),
                                make.default = TRUE) 
} else {
  stop("Registry already exists!")
}

meta <- read_tsv(here('data', paste0('meta_', tolower(disease),'.tsv')))

# make a function to run a single test on a single/multiple datasets
.f <- function(test, studies, disease, simba.loc){
  # load SIMBA
  library("tidyverse")
  devtools::load_all(simba.loc)
  message(test)
  message(disease)
  message(studies)
  if (disease=='CRC'){
    meta <- read_tsv(here('data', 'meta_crc.tsv'),col_types = cols())
    feat <- read.table(here('data', 'motus_crc_rel_meta.tsv'), sep='\t',
                           stringsAsFactors = FALSE, check.names = FALSE, 
                           quote = '', comment.char = '')
    feat <- as.matrix(feat)
    feat.count <- read.table(here('data', 'motus_crc_meta.tsv'), sep='\t',
                                 stringsAsFactors = FALSE, check.names = FALSE, 
                                 quote = '', comment.char = '')
    feat.count <- as.matrix(feat.count[rownames(feat),])
  } else if (disease=='CD'){
    meta <- read_tsv(here('data', 'meta_cd.tsv'),col_types = cols())
    feat <- read.table(here('data', 'motus_cd_rel_meta.tsv'), sep='\t',
                       stringsAsFactors = FALSE, check.names = FALSE, 
                       quote = '', comment.char = '')
    feat <- as.matrix(feat)
    feat.count <- read.table(here('data', 'motus_cd_meta.tsv'), sep='\t',
                             stringsAsFactors = FALSE, check.names = FALSE, 
                             quote = '', comment.char = '')
    feat.count <- as.matrix(feat.count[rownames(feat),])
  } else {
    stop("Unknown disease")
  }
  included.studies <- unique(meta$Study)
  studies <- as.logical(as.numeric(strsplit(studies, split = '')[[1]]))
  included.studies <- included.studies[studies]
  meta <- meta[meta$Study %in% included.studies,]
  feat.count <- feat.count[,meta$Sample_ID]
  
  # iterate over repetitions
  label <- recode(meta$Group, 'CTR'=-1, 'CRC'=1, 'CD'=1)
  names(label) <- meta$Sample_ID
  if (test %in% c('lm', 'wilcoxon')){
    feat.count <- SIMBA:::norm.data(feat.count, 'TSS.log', 1e-05)
  } else if (test=='limma'){
    feat.count <- SIMBA:::norm.data(feat.count, 'TSS.arcsin', 1e-05)
  }
  p.val <- SIMBA:::run.test(data=feat.count, label=label, test=test, conf=NULL)
  return(p.val)
}

# make combinations of studies for CRC
full.grid <- expand.grid(rep(list(0:1), length(unique(meta$Study))))
full.grid <- full.grid[rowSums(full.grid) > 1,]
combinations <- apply(full.grid, 1, paste, collapse='')

# add jobs
batchMap(.f, 'test'=tests, 'studies'=rep(combinations, each=length(tests)),
         more.args=list('disease'=disease, 'simba.loc'=simba.loc))

submitJobs()

# ##############################################################################
# to be executed manually once all jobs have finished
if (FALSE){
  message("collecting results")
  loadRegistry(paste0(temp.loc,  'test_results_registries/real_data_', disease))
  meta <- read_tsv(here('data', paste0('meta_', tolower(disease),'.tsv')),
                   col_types=cols())
  temp <- getJobPars() %>% 
    unnest_wider(job.pars)
  result.list <- list()
  pb <- progress::progress_bar$new(total=length(unique(temp$studies)))
  for (i in unique(temp$studies)){
    job.ids <- temp %>% filter(studies==i) %>% pull(job.id)
    p.val.mat <- reduceResults(cbind, job.ids)
    colnames(p.val.mat) <- temp %>% filter(studies==i) %>% pull(test)
    used.studies <- unique(meta$Study)[as.logical(
      as.numeric(strsplit(i, split='')[[1]]))] %>% 
      paste(collapse = '-')
    result.list[[used.studies]] <- p.val.mat
    pb$tick()
  }
  save(result.list, file=here('files', paste0('all_tests_', disease, '.Rdata')))
}
