if (!requireNamespace("renv", quietly = TRUE)){
  install.packages("renv")
}
renv::activate()

renv::restore()
