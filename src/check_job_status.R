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

reg <- here('test_results_registries', test.dir)

if (!dir.exists(reg)) {
	stop("Directory does not exist")
} 

loadRegistry(reg)

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

