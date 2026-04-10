library(tidyverse)

# CONFIG
HMM_DIR      <- "output/hmm/search"
ORF_DIR      <- "output/orfs"
MANIFEST_TSV <- "output/manifest_genomes.tsv"
OUT_DIR      <- "output/VMD-database/"

OUT_TABLES  <- file.path(OUT_DIR, "virus_informations")
OUT_MARKERS <- file.path(OUT_DIR, "markers")

dir.create(OUT_TABLES,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_MARKERS, recursive = TRUE, showWarnings = FALSE)

# MANIFEST
manifest <- read_tsv(MANIFEST_TSV, show_col_types = FALSE)

# read the tblout -> raw
tbl_files <- list.files(HMM_DIR, pattern = "\\.tbl$", full.names = TRUE)

read_tblout <- function(tbl) {
  lines <- readLines(tbl)
  lines <- lines[!startsWith(lines, "#")] 
  if (length(lines) == 0) return(tibble())

  fields <- strsplit(lines, "\\s+")

  tibble(
    virus_id = sub("__.*$", "", basename(tbl)),
    marker   = sub("^.*__", "", str_remove(basename(tbl),".tbl")),
    orf_name = vapply(fields, `[`, "", 1),  
    evalue   = as.numeric(vapply(fields, `[`, "", 5)), 
    score    = as.numeric(vapply(fields, `[`, "", 6)), 
    description = vapply(fields, function(x) paste(x[19:length(x)], collapse = " "), "") 
  )
}

tblout <- map_dfr(tbl_files, read_tblout)

# best hit by (virus_id, marker)
best_hits <- tblout %>%
  group_by(virus_id, marker) %>%
  arrange(evalue, desc(score)) %>%
  slice(1) %>%
  ungroup()

# Extract prodigal (start, end, strand)
parts <- strsplit(best_hits$description, "#", fixed = TRUE)

best_hits <- best_hits %>%
  mutate(
    prodigal_start = as.integer(trimws(sapply(parts, `[`, 2))),
    prodigal_end = as.integer(trimws(sapply(parts, `[`, 3)))
  ) %>%
  left_join(manifest, by = "virus_id")


# Export table

# General virus composition
# Pivot table to make it shorter 
markers <- sort(unique(best_hits$marker))
taxonomy_all <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
virus_compo_taxo <- best_hits %>%
  select(virus_id, Virus_names, marker, orf_name) %>%
  pivot_wider(
    names_from = marker,
    values_from = orf_name,
  ) %>%
  left_join(
    manifest %>%
     select(virus_id, ICTV_ID, all_of(taxonomy_all)),
    by = "virus_id"
  ) %>%
  select(virus_id, Virus_names, ICTV_ID, all_of(markers), all_of(taxonomy_all))

write_tsv(virus_compo_taxo, file.path(OUT_TABLES, "virus_compo_taxo.tsv"))

# Virus metadata
virus_metadata <- manifest %>%
  select(virus_id, ICTV_ID, Virus_names, Virus_names_abrv, Host_source,source, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  rename(Origin_source = source) %>%
  filter(virus_id %in% virus_compo_taxo$virus_id)

write_tsv(virus_metadata, file.path(OUT_TABLES, "virus_metadata.tsv"))

# Folder by marker
load_orfs <- function(virus_id) {
  faa_path <- file.path(ORF_DIR, str_c(virus_id, ".faa"))
  fna_path <- file.path(ORF_DIR, str_c(virus_id, ".fna"))

  aa_seq <- if (file.exists(faa_path)) Biostrings::readAAStringSet(faa_path) else Biostrings::AAStringSet()
  nt_seq <- if (file.exists(fna_path)) Biostrings::readDNAStringSet(fna_path) else Biostrings::DNAStringSet()

  names(aa_seq) <- sub(" .*", "", names(aa_seq))
  names(nt_seq) <- sub(" .*", "", names(nt_seq))

  list(aa = aa_seq, nt = nt_seq)
}

ORF_CACHE <- list()
get_orfs <- function(virus_id) {
  if (is.null(ORF_CACHE[[virus_id]])) {
    ORF_CACHE[[virus_id]] <- load_orfs(virus_id)
  }
  ORF_CACHE[[virus_id]]
}


# Remove trailing stop-codon indicator (*)
strip_terminal_star <- function(x) {
  x_chr <- as.character(x)
  x_chr2 <- sub("\\*$", "", x_chr)
  names(x_chr2) <- names(x)
  Biostrings::AAStringSet(x_chr2)
}



for (m in markers) {
  out_m <- file.path(OUT_MARKERS, m)
  dir.create(out_m, recursive = TRUE, showWarnings = FALSE)

  df_m <- best_hits %>%
    filter(marker == m) %>%
    arrange(evalue, desc(score))

  # Export manifest info 
  df_m2 <- df_m %>%
    select(virus_id,marker,orf_name,evalue,score,prodigal_start,prodigal_end,Species,Virus_names)
  write_tsv(df_m2, file.path(out_m, str_c(m,"_virus_orf_description.tsv")))

  aa_m <- Biostrings::AAStringSet()
  nt_m <- Biostrings::DNAStringSet()

  for (i in seq_len(nrow(df_m))) {
    r <- df_m[i, ]

    orfs <- get_orfs(r$virus_id)

    hdr <- str_c(
      r$virus_id,
      " ",
      df_m[i, ] %>%
        unite("Taxonomy", Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ";") %>%
        pull(Taxonomy)
      )

    aa_m <- c(aa_m, orfs$aa[r$orf_name])
    aa_m <- strip_terminal_star(aa_m) 
    names(aa_m)[length(aa_m)] <- hdr

    nt_m <- c(nt_m, orfs$nt[r$orf_name])
    names(nt_m)[length(nt_m)] <- hdr
   }

  # Write sequences FASTA files (AA and NT)
  Biostrings::writeXStringSet(aa_m, file.path(out_m, str_c(m,"_protein.fasta")), format = "fasta", width = 60)
  Biostrings::writeXStringSet(nt_m, file.path(out_m,  str_c(m,"_nucleotid.fasta")), format = "fasta", width = 60)

}

# Export database in specific tools format 
source("scripts/utils/export_format_database.R")
