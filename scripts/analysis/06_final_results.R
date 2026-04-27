library(tidyverse)

# CONFIG
HMM_DIR         <- "output/hmm/search"
ORF_DIR         <- "output/orfs"
REF_PROTEIN_DIR <- "output/references_protein"
MANIFEST_TSV    <- "output/config/manifest_genomes.tsv"
MAP_TSV         <- "output/config/marker_taxo_map.tsv"
OUT_DIR         <- "output/VMD-database"
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

# READ HMMSEARCH TABLES
tbl_files <- list.files(HMM_DIR, pattern = "\\.tbl$", full.names = TRUE)

read_tblout <- function(tbl) {
  lines <- readLines(tbl, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]

  if (length(lines) == 0) {
    return(tibble())
  }

  fields <- strsplit(trimws(lines), "\\s+")
  fname <- str_remove(basename(tbl), "\\.tbl$")

  name_parts <- strsplit(fname, "__", fixed = TRUE)[[1]]
  virus_id <- name_parts[1]
  marker_group_id <- paste(name_parts[-1], collapse = "__")

  tibble(
    virus_id = virus_id,
    marker_group_id = marker_group_id,
    orf_name = vapply(fields, `[`, "", 1),
    evalue = as.numeric(vapply(fields, `[`, "", 5)),
    score = as.numeric(vapply(fields, `[`, "", 6)),
    description = vapply(fields, function(x) paste(x[19:length(x)], collapse = " "), "")
  )
}

tblout <- map_dfr(tbl_files, read_tblout)

# Selection of ORF

default_threshold_score      <- 0.80
default_threshold_length     <- 0.80
default_ref_length_fraction  <- 0.75

# Extract prodigal coordinates
parts <- strsplit(tblout$description, "#", fixed = TRUE)

tblout <- tblout %>%
  mutate(
    prodigal_start = as.integer(trimws(sapply(parts, `[`, 2))),
    prodigal_end   = as.integer(trimws(sapply(parts, `[`, 3))),
    orf_length     = abs(prodigal_end - prodigal_start) + 1
  )

# Reference length thresholds
ref_length_thresholds <- list.files(
  REF_PROTEIN_DIR,
  pattern = "\\.fasta$",
  full.names = TRUE
) %>%
  map_dfr(function(path) {
    ref_seq <- Biostrings::readAAStringSet(path)
    ref_len <- as.integer(Biostrings::width(ref_seq))

    tibble(
      marker_group_id = tools::file_path_sans_ext(basename(path)),
      n_ref_sequences = length(ref_seq),
      min_ref_length_aa = min(ref_len),
      max_ref_length_aa = max(ref_len),
      min_detected_orf_length_aa = ceiling(min(ref_len) * default_ref_length_fraction)
    )
  })

write_tsv(
  ref_length_thresholds,
  file.path(OUT_LOGS, "ref_length_thresholds.tsv")
)

# ORF protein lengths
load_orf_lengths <- function(virus_id) {
  faa_path <- file.path(ORF_DIR, str_c(virus_id, ".faa"))

  if (!file.exists(faa_path)) {
    return(tibble())
  }

  aa_seq <- Biostrings::readAAStringSet(faa_path)

  tibble(
    virus_id = virus_id,
    orf_name = sub(" .*", "", names(aa_seq)),
    orf_length_aa = as.integer(Biostrings::width(aa_seq))
  )
}

orf_lengths <- map_dfr(unique(tblout$virus_id), load_orf_lengths)

# Remove short ORFs before ranking
tblout_all_detected <- tblout %>%
  left_join(
    orf_lengths,
    by = c("virus_id", "orf_name")
  ) %>%
  left_join(
    ref_length_thresholds %>%
      select(marker_group_id, min_ref_length_aa, min_detected_orf_length_aa),
    by = "marker_group_id"
  ) %>%
  mutate(
    pass_min_ref_length = orf_length_aa >= min_detected_orf_length_aa
  ) %>%
  filter(pass_min_ref_length)

tblout_removed_by_ref_length <- tblout %>%
  left_join(
    orf_lengths,
    by = c("virus_id", "orf_name")
  ) %>%
  left_join(
    ref_length_thresholds %>%
      select(marker_group_id, min_ref_length_aa, min_detected_orf_length_aa),
    by = "marker_group_id"
  ) %>%
  mutate(
    pass_min_ref_length = orf_length_aa >= min_detected_orf_length_aa
  ) %>%
  filter(!pass_min_ref_length)

write_tsv(
  tblout_removed_by_ref_length,
  file.path(OUT_LOGS, "orfs_removed_by_min_ref_length.tsv")
)

write_tsv(
  tblout_removed_by_ref_length %>%
    group_by(marker_group_id) %>%
    summarise(
      n_orfs_removed = n(),
      .groups = "drop"
    ),
  file.path(OUT_LOGS, "orfs_removed_by_min_ref_length_summary.tsv")
)

# Rank hits inside each virus x marker_group
tblout_ranked <- tblout_all_detected %>%
  group_by(virus_id, marker_group_id) %>%
  arrange(evalue, desc(score), .by_group = TRUE) %>%
  mutate(
    copy_rank = row_number(),
    best_score = first(score),
    best_orf_length = first(orf_length),
    score_ratio = score / best_score,
    length_ratio = orf_length / best_orf_length
  ) %>%
  ungroup()

# Keep best hit + probable extra copies
selected_orfs <- tblout_ranked %>%
  filter(
    copy_rank == 1 |
      (
        score_ratio >= default_threshold_score &
        length_ratio >= default_threshold_length
      )
  ) %>%
  left_join(marker_map, by = "marker_group_id") %>%
  left_join(manifest, by = "virus_id")

# EXPORT TABLES
marker_group_ids <- sort(unique(selected_orfs$marker_group_id))

virus_compo_taxo <- selected_orfs %>%
  group_by(virus_id, Virus_names, marker_group_id) %>%
  summarise(
    orf_names = str_c(sort(unique(orf_name)), collapse = ";"),
    n_copies = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = marker_group_id,
    values_from = c(orf_names, n_copies),
    names_glue = "{marker_group_id}_{.value}",
    values_fill = list(orf_names = NA_character_, n_copies = 0)
  ) %>%
  left_join(
    manifest %>%
      select(virus_id, ICTV_ID, all_of(taxonomy_all)),
    by = "virus_id"
  )

orf_cols <- str_c(marker_group_ids, "_orf_names")
copy_cols <- str_c(marker_group_ids, "_n_copies")

virus_compo_taxo <- virus_compo_taxo %>%
  select(
    virus_id, Virus_names, ICTV_ID,
    all_of(orf_cols),
    all_of(copy_cols),
    all_of(taxonomy_all)
  )

write_tsv(virus_compo_taxo, file.path(OUT_TABLES, "virus_compo_taxo.tsv"))

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

# Export all detected ORFs by marker
source("scripts/utils/export_all_detected_orfs_by_marker.R")

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
    arrange(virus_id, copy_rank, evalue, desc(score))

  df_m2 <- df_m %>%
    select(
      orf_name, virus_id, group_id, marker,
      copy_rank,
      evalue, score,
      best_score, score_ratio,
      orf_length, best_orf_length, length_ratio,
      prodigal_start, prodigal_end,
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
    file.path(out_m, str_c(mg, "_nucleotid.fasta")),
    format = "fasta",
    width = 60
  )
}

# Export database in specific tools format
source("scripts/utils/export_format_database.R")