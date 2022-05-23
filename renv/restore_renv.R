# ##############################################################################
#
## Restore the whole project with more memory
#
# ##############################################################################


renv::install('bioc::distinct')
renv::install('bioc::mixOmics')
renv::restore()
