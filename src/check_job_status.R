# ##############################################################################
#
# check status of job registry
#
# ##############################################################################

library("batchtools")
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("here"))

args <- commandArgs(trailingOnly = TRUE)
test.dir <- args[1]

reg <- paste0(temp.loc, 'test_results_registries/', test.dir)
message(reg)
if (!dir.exists(reg)) {
	stop("Directory does not exist")
} 

suppressMessages(loadRegistry(reg))

print(getStatus())

info <- getJobTable()
done.jobs <- findDone()

message("those tests are done already!")
print(info %>% filter(job.id %in% done.jobs$job.id) %>% pull(algorithm) %>% table())

done.jobs <- findNotDone()
message("those tests are not done year!")
print(info %>% filter(job.id %in% done.jobs$job.id) %>% pull(algorithm) %>% table())

done.jobs <- findErrors()
message("those tests errored!")
print(info %>% filter(job.id %in% done.jobs$job.id) %>% pull(algorithm) %>% table())

