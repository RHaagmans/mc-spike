renv_status <- renv::status()

if(!renv_status$synchronized){
  print("Renv not synchronized, trying restore")
  renv::restore()
}
