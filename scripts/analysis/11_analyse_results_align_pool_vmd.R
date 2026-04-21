# library(tidyverse)
# library(ggplot2)

# ncbi_dir <- "output/benchmark/pool_protein/pool_protein_raw"
# hits_dir <- "output/benchmark/alignment_vmd_pool/alignment"
# out_dir  <- "output/benchmark/alignment_vmd_pool/bench_results"
# missing_prot_path <- "output/benchmark/alignment_vmd_pool/missing_prot"
# hits_prot_path <- "output/benchmark/alignment_vmd_pool/hits_prot"

# dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# dir.create(missing_prot_path, recursive = TRUE, showWarnings = FALSE)
# dir.create(hits_prot_path, recursive = TRUE, showWarnings = FALSE)

# # read compo viro
# meta <- read_tsv("output/VMD-database/virus_informations/virus_compo_taxo.tsv")


# analyse_res <- function(marker,ncbi_dir,hits_dir,out_dir,missing_prot_path,hits_prot_path,meta){ 
#   prot <- Biostrings::readAAStringSet(str_c(ncbi_dir,"/", marker,".fasta"))

#   tab <- read.table(str_c(hits_dir,"/",marker,".m8"))
#   colnames(tab) <- c("query","ref","ident","alnlen","mismatch","gapopen","qstart","qend","sstart","send","evalue","bits") 

#   # Select best hit
#   tab_best <- tab |>
#     group_by(query)  |>
#     slice_min(evalue, with_ties = FALSE) 

#   # Write best match results
#   write_tsv(tab_best, str_c(out_dir,"/",marker,"_best_hits.tsv"))

#   # Let's check missing prot query 
#   prot_id <- sub(" .*", "", names(prot))
#   missing_prot_id <- setdiff(unique(prot_id), unique(tab$query))
#   prot_missing <- prot[prot_id %in% missing_prot_id]
  

#   # write in fasta missing prot 
#   if (length(prot_missing) > 0) {
#     Biostrings::writeXStringSet(prot_missing,str_c(missing_prot_path,"/",marker,".fasta"))
#   }

#   # Export hits prot 

#   prot_hits <- prot[prot_id %in% tab_best$query]
#   # write in fasta 
#   if (length(prot_hits) > 0) {
#     Biostrings::writeXStringSet(prot_hits,str_c(hits_prot_path,"/",marker,".fasta"))
#   }

#   # Visualise histogram of alignment
#   # Create histograms for each alignment quality metric
#   p1 <- ggplot(tab_best, aes(x = ident)) +
#     geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
#     labs(title = "Distribution of identity percentage (fident)",
#         x = "Identity percentage",
#         y = "Number of alignments") +
#     theme_minimal()

#   p2 <- ggplot(tab_best, aes(x = alnlen)) +
#     geom_histogram(bins = 30, fill = "darkgreen", color = "white", alpha = 0.8) +
#     labs(title = "Distribution of alignment lengths (alnlen)",
#         x = "Alignment length (aa)",
#         y = "Number of alignments") +
#     theme_minimal()

#   p3 <- ggplot(tab_best, aes(x = evalue)) +
#     geom_histogram(bins = 50, fill = "orange", color = "white", alpha = 0.8) +
#     scale_x_log10() +
#     labs(title = "Distribution of E-values (evalue)",
#         x = "E-value (log10 scale)",
#         y = "Number of alignments") +
#     theme_minimal()

#   p4 <- ggplot(tab_best, aes(x = bits)) +
#     geom_histogram(bins = 30, fill = "purple", color = "white", alpha = 0.8) +
#     labs(title = "Distribution of bit scores (bits)",
#         x = "Bit score",
#         y = "Number of alignments") +
#     theme_minimal()

#   # Arrange and display plots
#   final_plot <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)

#   # Save 
#   ggsave(
#     filename = str_c(out_dir,"/" ,marker,".svg"),
#     plot = final_plot,
#     width = 12,
#     height = 8,
#     units = "in",
#     dpi = 300)



#   # Analyse part of the database witch is not covered by the pool
#   # Family of the database ref with have the marker 
#   meta2 <- meta |>
#     filter(.data[[marker]] != "NA")
#   meta2$Family |> unique()

#   # Family of the reference database wich query blast 
#   meta3 <- meta2 |>
#     filter(virus_id %in% unique(tab$ref))
#   meta3$Family |> unique()
  


#   # Analyse on many query has best match with wich family 
  
#   family_count <- tab_best |>
#   left_join(meta, by = c("ref" = "virus_id")) |>
#   group_by(Family) |>
#   summarise(
#     nb_queries = n_distinct(query),  
#     .groups = "drop"
#   ) |>
#   arrange(desc(nb_queries)) 
#   print(family_count)


# }



# # Launch
# markers <- sub("\\.fasta$", "", list.files(ncbi_dir, pattern = "\\.fasta$", full.names = FALSE))

# lapply(markers, function(marker) analyse_res(marker,ncbi_dir,hits_dir,out_dir,missing_prot_path,hits_prot_path,meta))


























library(tidyverse)
library(ggplot2)

ncbi_dir <- "output/benchmark/pool_protein/pool_protein_raw"
hits_dir <- "output/benchmark/alignment_vmd_pool/alignment"
out_dir  <- "output/benchmark/alignment_vmd_pool/bench_results"
missing_prot_path <- "output/benchmark/alignment_vmd_pool/missing_prot"
hits_prot_path <- "output/benchmark/alignment_vmd_pool/hits_prot"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(missing_prot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(hits_prot_path, recursive = TRUE, showWarnings = FALSE)

meta <- read_tsv("output/VMD-database/virus_informations/virus_compo_taxo.tsv", show_col_types = FALSE)

analyse_res <- function(marker_group_id, ncbi_dir, hits_dir, out_dir, missing_prot_path, hits_prot_path, meta) {

  marker <- sub("^.*__", "", marker_group_id)

  pool_fasta <- file.path(ncbi_dir, str_c(marker, ".fasta"))
  hit_file <- file.path(hits_dir, str_c(marker_group_id, ".m8"))

  if (!file.exists(pool_fasta) || !file.exists(hit_file)) {
    return(invisible(NULL))
  }

  prot <- Biostrings::readAAStringSet(pool_fasta)

  tab <- read.table(hit_file)
  colnames(tab) <- c("query","ref","ident","alnlen","mismatch","gapopen","qstart","qend","sstart","send","evalue","bits")

  tab_best <- tab |>
    group_by(query) |>
    slice_min(evalue, with_ties = FALSE) |>
    ungroup()

  write_tsv(tab_best, file.path(out_dir, str_c(marker_group_id, "_best_hits.tsv")))

  prot_id <- sub(" .*", "", names(prot))
  missing_prot_id <- setdiff(unique(prot_id), unique(tab$query))
  prot_missing <- prot[prot_id %in% missing_prot_id]

  if (length(prot_missing) > 0) {
    Biostrings::writeXStringSet(prot_missing, file.path(missing_prot_path, str_c(marker_group_id, ".fasta")))
  }

  prot_hits <- prot[prot_id %in% tab_best$query]
  if (length(prot_hits) > 0) {
    Biostrings::writeXStringSet(prot_hits, file.path(hits_prot_path, str_c(marker_group_id, ".fasta")))
  }

  p1 <- ggplot(tab_best, aes(x = ident)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
    labs(title = "Distribution of identity percentage (fident)",
         x = "Identity percentage",
         y = "Number of alignments") +
    theme_minimal()

  p2 <- ggplot(tab_best, aes(x = alnlen)) +
    geom_histogram(bins = 30, fill = "darkgreen", color = "white", alpha = 0.8) +
    labs(title = "Distribution of alignment lengths (alnlen)",
         x = "Alignment length (aa)",
         y = "Number of alignments") +
    theme_minimal()

  p3 <- ggplot(tab_best, aes(x = evalue)) +
    geom_histogram(bins = 50, fill = "orange", color = "white", alpha = 0.8) +
    scale_x_log10() +
    labs(title = "Distribution of E-values (evalue)",
         x = "E-value (log10 scale)",
         y = "Number of alignments") +
    theme_minimal()

  p4 <- ggplot(tab_best, aes(x = bits)) +
    geom_histogram(bins = 30, fill = "purple", color = "white", alpha = 0.8) +
    labs(title = "Distribution of bit scores (bits)",
         x = "Bit score",
         y = "Number of alignments") +
    theme_minimal()

  final_plot <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)

  ggsave(
    filename = file.path(out_dir, str_c(marker_group_id, ".svg")),
    plot = final_plot,
    width = 12,
    height = 8,
    units = "in",
    dpi = 300
  )

  if (marker_group_id %in% colnames(meta)) {
    meta2 <- meta |>
      filter(.data[[marker_group_id]] != "NA")

    meta3 <- meta2 |>
      filter(virus_id %in% unique(tab$ref))

    family_count <- tab_best |>
      left_join(meta, by = c("ref" = "virus_id")) |>
      group_by(Family) |>
      summarise(
        nb_queries = n_distinct(query),
        .groups = "drop"
      ) |>
      arrange(desc(nb_queries))

    print(marker_group_id)
    print(unique(meta2$Family))
    print(unique(meta3$Family))
    print(family_count)
  }
}

marker_group_ids <- sub("\\.m8$", "", list.files(hits_dir, pattern = "\\.m8$", full.names = FALSE))

lapply(marker_group_ids, function(x) analyse_res(x,ncbi_dir,hits_dir,out_dir,missing_prot_path,hits_prot_path,meta))