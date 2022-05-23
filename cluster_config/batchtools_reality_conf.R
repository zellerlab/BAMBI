cluster.functions=makeClusterFunctionsSlurm(
  template=here('cluster_config', 'slurm_reality_check.tmpl'))
work.dir=here()
packages = c('renv', 'here')
source = c(here::here('src','load_pkgs.R'))
