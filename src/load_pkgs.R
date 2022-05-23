# #######################################################
#
## script sourced by batchtools registry on slave nodes
#
# #######################################################

# load the project-local library
renv::load(here::here())

# install local version of SIMBA
devtools::load_all(simba.loc)
