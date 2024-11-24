args<- commandArgs(trailingOnly = TRUE)

if(is.na(args[1])){
  # If no output path is supplied, run from the script location
  folder_dir<-"."
  output_dir<-".."
}else{
  folder_dir<-args[1]
  output_dir<-args[1]
}

mkdir <- function(subdir) {
  dir.create(
    file.path(folder_dir, subdir),
    showWarnings = FALSE,
    recursive = TRUE
  )
}

proj_report_path <- function(report_path){
  file.path(output_dir, "reports", report_path)
}

render_report <- function(rmd, report_filename){
  rmarkdown::render(
    rmd,
    output_file = proj_report_path(report_filename),
    params = list(output_dir = output_dir)
  )
  
}

mkdir("reports")
mkdir("figures")

mkdir("figures/read_and_contig_stats")
render_report('R/read_and_assembly_stats.Rmd','read_and_assembly_stats.html')

mkdir("figures/mc_bias")
render_report('R/mc_readmap_analysis.Rmd', 'mc_readmap_analysis.html')
render_report('R/mc_sequencing_bias.Rmd','mc_sequencing_bias.html')

mkdir("figures/mc_assembly")
render_report('R/mc_assembly_analysis.Rmd','mc_assembly_analysis.html')

mkdir("figures/virome_variation")
render_report(
  'R/virome_relative_abundance_variation.Rmd', 
  'virome_relative_abundance_variation.html'
)

mkdir("figures/virome_taxonomy")
render_report(
  'R/virome_taxonomy_diversity.Rmd',
  'virome_taxonomy_diversity.html'
)
