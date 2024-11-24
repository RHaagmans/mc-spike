library(tidyverse)

args<- commandArgs(trailingOnly = TRUE)

if(!is.na(args[1])){
  setwd(args[1])
}


read_stats <- read_tsv("data/read_qc/readstats.txt") %>%
  separate(file, into = c("base_sample", "spike", "rt_method", "read_pair")) %>%
  group_by(filter_step, base_sample, spike, rt_method) %>%
  summarise(
    n_reads = sum(num_seqs),
    total_bp = sum(sum_len),
    mean_length = mean(avg_len),
    mean_gc = mean(`GC(%)`)
  ) %>%
  mutate(
    library_name = paste0(base_sample, "-", spike, "-", rt_method),
    sample_name = paste0(base_sample, "-", spike)
  ) %>%
  relocate(library_name, sample_name, .after = filter_step)

write_tsv(read_stats, file = "data/read_qc/processed_read_stats.txt")

process_quast_data <- function(quast_file){
  read_tsv(quast_file) %>%
    separate(Assembly, into = c("assembly", "stage"), sep = "_contigs") %>%
    mutate(
      stage = na_if(stage, ""),
      stage = replace_na(stage, "assembly"),
      stage = str_remove(stage, "^\\.")
    ) %>%
    rename(
      n_contigs = `# contigs`,
      total_length = `Total length`,
      n_contigs_10k = `# contigs (>= 10000 bp)`
    ) %>%
    select(assembly, stage, n_contigs, n_contigs_10k, total_length, N50, L50)
}

ass_stats <- process_quast_data(
    "data/quast/assembly/quast_results/latest/combined_reference/transposed_report.tsv"
  )
write_tsv(
  ass_stats, 
  file = "data/assemblies/quast-assembly-stats_assembly.txt"
  )

coas_stats <- process_quast_data(
    "data/quast/coassembly/quast_results/latest/combined_reference/transposed_report.tsv"
  )
write_tsv(
  coas_stats, 
  file = "data/assemblies/quast-assembly-stats_coassembly.txt"
  )
