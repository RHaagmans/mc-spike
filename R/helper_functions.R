library(tidyverse)
library(phyloseq)
source("plot_settings.R")

add_assembly_meta <- function(df, sample_sheet){
  df %>% 
    left_join(
      sample_sheet %>% 
        select(-sample_name, -coassembly, -replicate) %>%  
        distinct(),
      by = join_by(assembly)
    ) 
}

add_coassembly_meta <- function(df, sample_sheet){
  df %>% 
    left_join(
      sample_sheet %>% 
        select(-sample_name, -assembly, -spike, -replicate) %>%  
        distinct(),
      by = join_by(assembly == coassembly)
    ) 
}

load_sample_sheet <- function(sample_sheet_file){
  read_csv(
    sample_sheet_file,
    col_types = cols(
      sample_name = col_factor(levels = full_sample_order),
      base_sample = col_factor(levels = base_sample_order),
      spike = col_factor(levels = spike_order),
      rt_method = col_factor(levels = method_order),
      replicate = col_factor(),
      assembly = col_factor(),
      coassembly = col_factor(),
      fw_readfile = col_character(),
      rv_readfile = col_character(),
    )
  )
}

load_quast_stats <- function(quast_stats_file){
  read_tsv(
    quast_stats_file,
    col_types = cols(
      assembly = col_character(),
      stage = col_factor(levels = c("assembly", "derep", "vir")),
      n_contigs = col_integer(),
      n_contigs_10k = col_integer(),
      total_length = col_integer(),
      N50 = col_integer(),
      L50 = col_integer()
    )) %>% 
    arrange(assembly)
}

load_contig_stats <- function(contig_stats_file){
  read_tsv(
    contig_stats_file,
    col_types = cols(
      assembly_name = col_factor(),
      contig = col_character(), 
      length = col_integer(),
      gc = col_double()
    ) 
  ) %>%
    rename(assembly = assembly_name) %>% 
    group_by(assembly) %>% 
    arrange(desc(length)) %>% 
    mutate(
      contig_idx = row_number(), 
      cum_length = cumsum(length)
    )
}

process_idxstats <- function(idxstats) {
  return(
    left_join(
      idxstats %>% filter(ref != "*"),
      left_join(
        idxstats %>%
          filter(ref == "*") %>%
          mutate(tot_unmapped = unmapped_reads) %>%
          select(assembly, sample, tot_unmapped),
        idxstats %>%
          group_by(assembly, sample, rt_method) %>%
          summarise(tot_mapped = sum(mapped_reads)) %>%
          select(assembly, sample, tot_mapped),
        by = c('assembly', 'sample')
      ),
      by =  c('assembly', 'sample')
    ) %>%
      group_by(assembly, sample) %>% 
      mutate(
        tot_reads = tot_unmapped + tot_mapped,
        frac_reads = mapped_reads / tot_reads,
        rpkm = (mapped_reads / (length/10^3)) / (tot_reads/10^6),
        relab = rpkm/sum(rpkm)
      ) %>% 
      ungroup()
  )
}

load_idxstats <- function(path, set_assembly = FALSE) {
  idxstats <- read_tsv(
    path,
    col_types = cols(
      assembly = col_factor(),
      sample = col_factor(levels = full_sample_order),
      base_sample = col_factor(levels = base_sample_order),
      spike = col_factor(levels = spike_order),
      rt_method = col_factor(levels = method_order),
      replicate = col_factor(),
      ref = col_factor(),
      length = col_number(),
      mapped_reads = col_number(),
      unmapped_reads = col_number()
    )
  )
  
  if(is.character(set_assembly)) {
    idxstats <- idxstats %>% mutate(assembly = set_assembly)
  } else if (!"assembly" %in% names(idxstats)) {
    idxstats <- idxstats %>% mutate(assembly = "target_sequence")
  }
  
  return(process_idxstats(idxstats))
}

summarise_idxstats <- function(idxstats) {
  return(
    idxstats %>%
      group_by(assembly, sample, base_sample, spike, rt_method) %>%
      summarise(
        total = mean(tot_reads),
        mapped = mean(tot_mapped),
        .groups = 'keep'
      ) %>%
      mutate(percent = mapped / total,
             per_10k = mapped / (total / 10 ^ 4))
  )
}

sum_segment_readcounts <- function(idxstats) {
  return(
    idxstats %>%
      group_by(sample, base_sample, spike, rt_method, virus) %>%
      summarise(
        ref_length = sum(length),
        mapped_reads = sum(mapped_reads),
        tot_mapped = mean(tot_mapped),
        tot_reads = mean(tot_reads)
      )
  )
}

calc_expected_levels <- function(idxstats_by_virus,
                                 virus_quantities,
                                 correct_for_strands = TRUE) {
    idxstats_by_virus %>%
      left_join(virus_quantities,
                by = "virus") %>%
      group_by(sample, rt_method) %>%
      mutate(
        strand_correction = correct_for_strands,
        mc_expected_fraction = if_else(
          strand_correction,
          (ref_length * particles_added * strands) / sum(ref_length * particles_added * strands),
          (ref_length * particles_added) / sum(ref_length * particles_added)
        )
      ) %>%
      select(-strand_correction) %>%
      ungroup()
}

filter_faecal_virus_sequences <- function(idxstats, genomad_results) {
  idxstats %>%
    semi_join(genomad_results, by = c("assembly", "contig")) %>%
    # Recalculcate relative abundance only based on viral contigs
    group_by(assembly, sample) %>%
    mutate(rpkm = (mapped_reads / (length / 10 ^ 3)) / (tot_reads / 10 ^ 6),
           relab = rpkm / sum(rpkm)) %>%
    ungroup()
}

transform_idxstats_to_matrix <- function(idxstats, data_column){
  return(
    as.matrix(
      idxstats %>%    
        # R replaces "-" with "." when making df's and matrices...
        mutate(sample = str_replace_all(sample, "-", "_"),
               unique_contig = paste0(assembly, "_", contig)) %>%
        select(unique_contig, sample, all_of(data_column)) %>%
        pivot_wider(
          names_from = unique_contig, 
          values_from = all_of(data_column)) %>%
        remove_rownames() %>% 
        column_to_rownames(var="sample") %>%
        mutate(across(everything(), ~replace_na(.x, 0)))
    )
  )
}

phyloseq_to_tibble <- function(ps, sample_meta){
  otus <- as_tibble(
    as.data.frame(t(otu_table(ps, taxa_are_rows = FALSE))), 
    rownames = "unique_contig") %>%
    pivot_longer(
      cols = -unique_contig, names_to = "sample", values_to = "relab"
    ) %>%
    left_join(
      sample_meta,
      by = "sample"
    )
  tax <- as_tibble(as.data.frame(tax_table(ps)), rownames = "unique_contig")
  data <- left_join(
    otus,
    tax,
    by="unique_contig"
  )
  
  return(data)
}

load_coverage <- function(coverage_file){
  read_tsv(coverage_file,
           col_types = cols(
             assembly = col_factor(),
             sample = col_factor(levels = full_sample_order),
             base_sample = col_factor(levels = base_sample_order),
             spike = col_factor(levels = spike_order),
             rt_method = col_factor(levels = method_order),
             replicate = col_factor(),
             `#rname` = col_factor(),
             startpos = col_integer(),
             endpos = col_integer(),
             numreads = col_integer(),
             covbases = col_integer(),
             coverage = col_double(),
             meandepth = col_double(),
             meanbaseq = col_double(),
             meanmapq = col_double()
           )) %>% 
    rename(contig = `#rname`)
}

load_depth_table <- function(depth_stats_file){
  depth_table_data <- read_tsv(
    depth_stats_file,
    col_types = cols(
      sample=col_factor(levels = full_sample_order),
      base_sample=col_factor( levels = base_sample_order),
      spike=col_factor(levels = spike_order),
      rt_method=col_factor(levels = method_order),
      ref=col_factor(levels = segment_order),
      bp=col_number(),
      depth=col_number()
    ),
    show_col_types = FALSE) 
  return(
    depth_table_data %>%
      group_by(sample, rt_method, ref) %>% 
      mutate(depth_norm = depth/sum(depth)) %>%
      ungroup()
  )
}

calculate_coverage <- function(depth_table){
  return(
    depth_table %>% 
      group_by(sample, base_sample, spike, rt_method, ref) %>% 
      summarise(
        total_bp = n(), 
        covered = sum(depth>0), 
        coverage = covered/total_bp
      ) %>% 
      ungroup()
  )
  
  
}

load_sequence_stats <- function(sequence_stats_file) {
  return(read_tsv(
    sequence_stats_file,
    col_names = c("ref", "ref_length", "ref_gc"),
    col_types = cols(
      ref = col_factor(),
      ref_length = col_integer(),
      ref_gc = col_number()
    ),
    show_col_types = FALSE
  ))
}

load_blast_ani_results <- function(blast_results_file){
  read_tsv(
    blast_results_file,
    col_types = cols(
      assembly_name = col_factor(),
      contig = col_factor(),
      ref = col_factor(),
      num_alns = col_number(),
      ani = col_number(),
      qcov = col_number(),
      tcov = col_number()
    )
  ) %>%
    rename(assembly = assembly_name)
  
  
}


load_genomad <- function(path) {
  return(
    read_tsv(
      path,
      col_types = cols(
        assembly_name = col_factor(),
        seq_name = col_factor(),
        length = col_integer(),
        topology = col_factor(),
        coordinates = col_character(),
        n_genes = col_number(),
        genetic_code = col_integer(),
        virus_score = col_number(),
        fdr = col_number(),
        n_hallmarks = col_number(),
        marker_enrichment = col_number(),
        taxonomy = col_character()
      )
    ) %>%
      rename(assembly = assembly_name,
             contig = seq_name)
  )
}

load_checkv_data <- function(checkv_file){
  read_tsv(
    checkv_file,
    col_types = cols(
      assembly_name = col_factor(),
      contig_id = col_character(),
      contig_length = col_integer(),
      provirus = col_factor(),
      proviral_length = col_integer(),
      gene_count = col_integer(),
      viral_genes = col_integer(),
      host_genes = col_integer(),
      checkv_quality = col_factor(levels = names(checkv_quality_palette)),
      miuvig_quality = col_factor(),
      completeness = col_double(),
      completeness_method = col_factor(),
      contamination = col_double(),
      kmer_freq = col_double(),
      warnings = col_character()
    )) %>% 
    rename(
      assembly = assembly_name,
      contig = contig_id
      ) %>% 
    mutate(
      checkv_quality_num = case_match(
        checkv_quality,
        "Complete"~ 4,
        "High-quality" ~ 3,
        "Medium-quality" ~ 2,
        "Low-quality" ~ 1,
        "Not-determined" ~ 0
      )
    )
}
