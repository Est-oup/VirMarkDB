library(tidyverse)

# CONFIG
HMM_DIR         <- "output/hmm/search"
ORF_DIR         <- "output/orfs"
REF_PROTEIN_DIR <- "output/references_protein"
MANIFEST_TSV    <- "output/config/manifest_genomes.tsv"
MAP_TSV         <- "output/config/marker_taxo_map.tsv"
OUT_DIR         <- "output/VirMarkDB"
OUT_LOGS        <- "output/hmm/logs"

OUT_TABLES      <- file.path(OUT_DIR, "virus_informations")
OUT_MARKERS     <- file.path(OUT_DIR, "markers")

dir.create(OUT_TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_MARKERS, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_LOGS, recursive = TRUE, showWarnings = FALSE)

taxonomy_all <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# INPUTS
manifest <- read_tsv(MANIFEST_TSV, show_col_types = FALSE)

marker_map <- read_tsv(MAP_TSV, show_col_types = FALSE) %>%
  mutate(
    marker_group_id = if_else(
      is.na(marker_group_id) | marker_group_id == "",
      str_c(group_id, "__", marker),
      marker_group_id
    )
  ) %>%
  distinct(marker_group_id, group_id, marker)

# READ HMMSEARCH DOMTBLOUT TABLES
tbl_files <- list.files(HMM_DIR, pattern = "\\.tbl$", full.names = TRUE)

domtblout_cols <- c("orf_name","target_accession","target_length_aa","marker_group_id","hmm_accession","hmm_length_aa",
  "evalue","score","bias","domain_number","domain_total","c_evalue","i_evalue","domain_score","domain_bias","hmm_from",
  "hmm_to","ali_from","ali_to","env_from","env_to","acc"
)

read_domtblout <- function(tbl) {
  lines <- readLines(tbl, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]

  if (length(lines) == 0) {return(tibble())}

  hmmer_raw <- readr::read_table(
    str_c(str_c(lines, collapse = "\n"),"\n"),
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE
  )

  hmmer <- hmmer_raw %>%
    select(1:22)
  names(hmmer) <- domtblout_cols

  hmmer$description <- apply(hmmer_raw[, 23:ncol(hmmer_raw), drop = FALSE], 1,function(x) str_c(x, collapse = " "))

  hmmer <- hmmer %>%
    mutate(
      virus_id = str_remove(orf_name, "_[0-9]+$"),
      orf_index = str_extract(orf_name, "[0-9]+$"),
      .before = 1
    ) %>%
    mutate(
      across(c(target_length_aa,hmm_length_aa,domain_number,domain_total,hmm_from,hmm_to,ali_from,ali_to,env_from,env_to,orf_index),as.integer),
      across(c(evalue,score,bias,c_evalue,i_evalue,domain_score,domain_bias,acc),as.numeric),
      virus_id_cut = str_remove(virus_id,"\\.[0-9]+$")
    ) %>%
    select(virus_id,virus_id_cut, orf_name, target_length_aa, marker_group_id, hmm_length_aa, evalue, score, description) %>%
    distinct()
}

tblout <- map_dfr(tbl_files, read_domtblout)

# Keep best hit and credible additional copies
score_threshold_multicopy <- 0.8
length_threshold_multicopy <- 0.8

# Rank hits inside each virus x marker_group
tblout_ranked <- tblout %>%
  group_by(virus_id, marker_group_id) %>%
  arrange(evalue, desc(score), .by_group = TRUE) %>%
  mutate(
    copy_rank = row_number(),
    best_score = first(score),
    best_orf_length = first(target_length_aa),
    score_ratio = score / best_score,
    length_ratio = target_length_aa / best_orf_length
  ) %>%
  rename(virus_id_tblout = virus_id) %>%
  ungroup()

# Manage virus id in manifest
manifest <- manifest %>% 
  mutate(
    virus_id_cut = str_remove(virus_id, "_partial"),
  ) 

# Keep best hit + probable extra copies
selected_orfs <- tblout_ranked %>%
  filter(
    copy_rank == 1 |
      (
        score_ratio >= score_threshold_multicopy &
        length_ratio >= length_threshold_multicopy
      )
  ) %>%
  left_join(marker_map, by = "marker_group_id") %>%
  left_join(manifest, by = "virus_id_cut")

# EXPORT TABLES
marker_group_ids <- sort(unique(selected_orfs$marker_group_id))


# Write table with virus composition across markers + taxonomy
virus_compo_taxo <- selected_orfs %>%
  group_by(virus_id, marker_group_id) %>%
  summarise(
    orf_names = str_c(sort(unique(orf_name)), collapse = ";"),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = marker_group_id,
    values_from = orf_names,
    names_glue = "{marker_group_id}_{.value}",
    values_fill = NA_character_
  ) %>%
  left_join(
    manifest %>%
      select(virus_id, ICTV_ID, Virus_names, all_of(taxonomy_all)),
    by = "virus_id"
  )

orf_cols <- str_c(marker_group_ids, "_orf_names")

virus_compo_taxo <- virus_compo_taxo %>%
  select(
    virus_id, ICTV_ID, Virus_names,
    all_of(orf_cols),
    all_of(taxonomy_all)
  )

write_tsv(virus_compo_taxo, file.path(OUT_TABLES, "virus_compo_taxo.tsv"))


# Write table with virus metadata (specificity, sources, etc)

virus_metadata <- manifest %>%
  select(
    virus_id, ICTV_ID, Virus_names, Virus_names_abrv, Host_source, Source,
    Kingdom, Phylum, Class, Order, Family, Genus, Species
  ) %>%
  rename(Origin_source = Source) %>%
  filter(virus_id %in% virus_compo_taxo$virus_id)

write_tsv(virus_metadata, file.path(OUT_TABLES, "virus_metadata.tsv"))

# Calculate logs
source("scripts/utils/hmm_search_logs.R")

# LOAD ORFS
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

# Remove terminal star
strip_terminal_star <- function(x) {
  x_chr <- as.character(x)
  x_chr2 <- sub("\\*$", "", x_chr)
  names(x_chr2) <- names(x)
  Biostrings::AAStringSet(x_chr2)
}

# EXPORT FASTA BY GROUP / MARKER
targets <- selected_orfs %>%
  distinct(group_id, marker, marker_group_id) %>%
  arrange(group_id, marker)

for (i in seq_len(nrow(targets))) {
  g  <- targets$group_id[i]
  m  <- targets$marker[i]
  mg <- targets$marker_group_id[i]

  out_m <- file.path(OUT_MARKERS, g, m)
  dir.create(out_m, recursive = TRUE, showWarnings = FALSE)

  df_m <- selected_orfs %>%
    filter(marker_group_id == mg) %>%
    arrange(virus_id, copy_rank, evalue, desc(score)) %>%
    mutate(
      description = str_remove(description, "^# ")
    ) %>%
    separate(
      description,
      into = c("position_start", "position_end", "prodigal_strand", "prodigal_info"),
      sep = " # ",
      extra = "merge",
      convert = TRUE
    ) 

  df_m2 <- df_m %>%
    select(
      orf_name, virus_id, group_id, marker,
      copy_rank,
      evalue, score,
      best_score, score_ratio,
      target_length_aa, best_orf_length, length_ratio,
      position_start, position_end,
      Species, Virus_names
    )

  write_tsv(
    df_m2,
    file.path(out_m, str_c(mg, "_virus_orf_description.tsv"))
  )

  aa_m <- Biostrings::AAStringSet()
  nt_m <- Biostrings::DNAStringSet()

  for (j in seq_len(nrow(df_m))) {
    r <- df_m[j, ]

    orfs <- get_orfs(r$virus_id)

    hdr <- str_c(
      r$orf_name,
      " ",
      m,
      " ",
      df_m[j, ] %>%
        unite("Taxonomy", Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ";") %>%
        pull(Taxonomy)
    )

    aa_m <- c(aa_m, orfs$aa[r$orf_name])
    aa_m <- strip_terminal_star(aa_m)
    names(aa_m)[length(aa_m)] <- hdr

    nt_m <- c(nt_m, orfs$nt[r$orf_name])
    names(nt_m)[length(nt_m)] <- hdr
  }

  Biostrings::writeXStringSet(
    aa_m,
    file.path(out_m, str_c(mg, "_protein.fasta")),
    format = "fasta",
    width = 60
  )

  Biostrings::writeXStringSet(
    nt_m,
    file.path(out_m, str_c(mg, "_nucleotide.fasta")),
    format = "fasta",
    width = 60
  )
}

# Export database in specific tools format
source("scripts/utils/export_format_database.R")