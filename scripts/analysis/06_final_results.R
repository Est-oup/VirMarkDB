library(tidyverse)

# CONFIG
HMM_DIR      <- "output/hmm/search"
ORF_DIR      <- "output/orfs"
MANIFEST_TSV <- "output/config/manifest_genomes.tsv"
MAP_TSV      <- "output/config/marker_taxo_map.tsv"
OUT_DIR      <- "output/VMD-database"

OUT_TABLES  <- file.path(OUT_DIR, "virus_informations")
OUT_MARKERS <- file.path(OUT_DIR, "markers")

dir.create(OUT_TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_MARKERS, recursive = TRUE, showWarnings = FALSE)

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

# BEST HIT BY GENOME x MARKER_GROUP
best_hits <- tblout %>%
  group_by(virus_id, marker_group_id) %>%
  arrange(evalue, desc(score)) %>%
  slice(1) %>%
  ungroup()

# EXTRACT PRODIGAL START / END
parts <- strsplit(best_hits$description, "#", fixed = TRUE)

best_hits <- best_hits %>%
  mutate(
    prodigal_start = as.integer(trimws(sapply(parts, `[`, 2))),
    prodigal_end   = as.integer(trimws(sapply(parts, `[`, 3)))
  ) %>%
  left_join(marker_map, by = "marker_group_id") %>%
  left_join(manifest, by = "virus_id")

# EXPORT TABLES
marker_group_ids <- sort(unique(best_hits$marker_group_id))

virus_compo_taxo <- best_hits %>%
  select(virus_id, Virus_names, marker_group_id, orf_name) %>%
  pivot_wider(
    names_from = marker_group_id,
    values_from = orf_name
  ) %>%
  left_join(
    manifest %>%
      select(virus_id, ICTV_ID, all_of(taxonomy_all)),
    by = "virus_id"
  ) %>%
  select(virus_id, Virus_names, ICTV_ID, all_of(marker_group_ids), all_of(taxonomy_all))

write_tsv(virus_compo_taxo, file.path(OUT_TABLES, "virus_compo_taxo.tsv"))

virus_metadata <- manifest %>%
  select(
    virus_id, ICTV_ID, Virus_names, Virus_names_abrv, Host_source, Source,
    Kingdom, Phylum, Class, Order, Family, Genus, Species
  ) %>%
  rename(Origin_source = Source) %>%
  filter(virus_id %in% virus_compo_taxo$virus_id)

write_tsv(virus_metadata, file.path(OUT_TABLES, "virus_metadata.tsv"))
write_tsv(best_hits, file.path(OUT_TABLES, "best_hits.tsv"))

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

# REMOVE TERMINAL STAR
strip_terminal_star <- function(x) {
  x_chr <- as.character(x)
  x_chr2 <- sub("\\*$", "", x_chr)
  names(x_chr2) <- names(x)
  Biostrings::AAStringSet(x_chr2)
}

# EXPORT FASTA BY GROUP / MARKER
targets <- best_hits %>%
  distinct(group_id, marker, marker_group_id) %>%
  arrange(group_id, marker)

for (i in seq_len(nrow(targets))) {
  g  <- targets$group_id[i]
  m  <- targets$marker[i]
  mg <- targets$marker_group_id[i]

  out_m <- file.path(OUT_MARKERS, g, m)
  dir.create(out_m, recursive = TRUE, showWarnings = FALSE)

  df_m <- best_hits %>%
    filter(marker_group_id == mg) %>%
    arrange(evalue, desc(score))

  df_m2 <- df_m %>%
    select(
      virus_id, group_id, marker, marker_group_id,
      orf_name, evalue, score, prodigal_start, prodigal_end,
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

# EXPORT DATABASE IN SPECIFIC TOOLS FORMAT
source("scripts/utils/export_format_database.R")

