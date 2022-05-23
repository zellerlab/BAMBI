cluster.functions = makeClusterFunctionsSlurm(template = here('cluster_config', 'slurm_test.tmpl'))
# cluster.functions = makeClusterFunctionsMulticore(ncpus = future::availableCores()-1)
# cluster.functions = makeClusterFunctionsSGE(template = here('cluster_config', 'sge_test.tmpl'))
# cluster.functions = makeClusterFunctionsInteractive()

work.dir = here()
packages = c('renv', 'here') # required to load most up-to-date SIMBA on slaves

# bootstrap the renv project library on the slaves
source = c(here::here('src','load_pkgs.R'))

# enable memory measurement
default.resources = list(measure.memory = TRUE)