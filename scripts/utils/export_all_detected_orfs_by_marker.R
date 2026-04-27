library(tidyverse)

OUT_MARKERS_ALL <- file.path(OUT_LOGS, "markers_all_detected")

dir.create(OUT_MARKERS_ALL, recursive = TRUE, showWarnings = FALSE)

# Add metadata to all detected ORFs
all_detected_orfs <- tblout_all_detected %>%
  left_join(
    marker_map %>%
      select(marker_group_id, group_id, marker) %>%
      distinct(),
    by = "marker_group_id"
  ) %>%
  left_join(
    manifest %>%
      select(virus_id, Virus_names, any_of(taxonomy_all)),
    by = "virus_id"
  )

# Export all detected ORFs by marker
targets_all <- all_detected_orfs %>%
  distinct(group_id, marker, marker_group_id) %>%
  arrange(group_id, marker)

for (i in seq_len(nrow(targets_all))) {
  g  <- targets_all$group_id[i]
  m  <- targets_all$marker[i]
  mg <- targets_all$marker_group_id[i]

  out_m <- file.path(OUT_MARKERS_ALL, g, m)
  dir.create(out_m, recursive = TRUE, showWarnings = FALSE)

  df_m <- all_detected_orfs %>%
    filter(marker_group_id == mg) %>%
    arrange(virus_id, evalue, desc(score))

  df_m2 <- df_m %>%
    select(
      orf_name, virus_id, group_id, marker, marker_group_id,
      evalue, score,
      prodigal_start, prodigal_end,
      orf_length, orf_length_aa,
      min_ref_length_aa, min_detected_orf_length_aa,
      Species, Virus_names
    )

  write_tsv(
    df_m2,
    file.path(out_m, str_c(mg, "_all_detected_orf_description.tsv"))
  )

  aa_m <- Biostrings::AAStringSet()
  nt_m <- Biostrings::DNAStringSet()

  for (j in seq_len(nrow(df_m))) {
    r <- df_m[j, ]

    orfs <- get_orfs(r$virus_id)

    if (!r$orf_name %in% names(orfs$aa)) {
      next
    }

    if (!r$orf_name %in% names(orfs$nt)) {
      next
    }

    hdr <- str_c(
      r$orf_name,
      " ",
      df_m[j, ] %>%
        unite("Taxonomy", Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ";") %>%
        pull(Taxonomy)
    )

    aa_one <- orfs$aa[r$orf_name]
    aa_one <- strip_terminal_star(aa_one)
    names(aa_one) <- hdr
    aa_m <- c(aa_m, aa_one)

    nt_one <- orfs$nt[r$orf_name]
    names(nt_one) <- hdr
    nt_m <- c(nt_m, nt_one)
  }

  Biostrings::writeXStringSet(
    aa_m,
    file.path(out_m, str_c(mg, "_all_detected_protein.fasta")),
    format = "fasta",
    width = 60
  )

  Biostrings::writeXStringSet(
    nt_m,
    file.path(out_m, str_c(mg, "_all_detected_nucleotid.fasta")),
    format = "fasta",
    width = 60
  )
}