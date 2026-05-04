library(tidyverse)
library(fs)

outdir <- "output/VirMarkDB/export_format"
virus_compo_taxo_path <- "output/VirMarkDB/virus_informations/virus_compo_taxo.tsv"
map_tsv <- "output/config/marker_taxo_map.tsv"

taxonomy_all <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

dir_create(file.path(outdir, "vsearch"))
dir_create(file.path(outdir, "dada2"))

# Read input tables
virus_compo_taxo <- read_tsv(virus_compo_taxo_path, show_col_types = FALSE)

marker_map <- read_tsv(map_tsv, show_col_types = FALSE) %>%
  mutate(
    marker_group_id = if_else(
      is.na(marker_group_id) | marker_group_id == "",
      str_c(group_id, "_", marker),
      marker_group_id
    )
  ) %>%
  distinct(group_id, marker, marker_group_id)

# Build taxonomy file for vsearch
orf_cols <- grep("_orf_names$", names(virus_compo_taxo), value = TRUE)

taxo_vsearch <- virus_compo_taxo %>%
  select(any_of(c("virus_id", taxonomy_all, orf_cols))) %>%
  pivot_longer(
    cols = all_of(orf_cols),
    names_to = "marker_group_col",
    values_to = "orf_names",
    values_drop_na = TRUE
  ) %>%
  separate_rows(orf_names, sep = ";") %>%
  filter(!is.na(orf_names), orf_names != "") %>%
  mutate(reference_id = orf_names) %>%
  unite(
    taxonomy,
    any_of(taxonomy_all),
    sep = ";",
    remove = TRUE,
    na.rm = FALSE
  ) %>%
  select(reference_id, taxonomy) %>%
  distinct()

write_tsv(taxo_vsearch, file.path(outdir, "vsearch", "taxonomy.tsv"))

# Extract taxonomy from FASTA headers: "virus_id Kingdom;Phylum;...;Species"
get_taxo_from_names <- function(x) {
  parts <- str_split_fixed(names(x), " ", 3)
  taxo <- parts[, 3]
  taxo[taxo == ""] <- "NA;NA;NA;NA;NA;NA;NA"
  taxo
}

for (i in seq_len(nrow(marker_map))) {
  g  <- marker_map$group_id[i]
  m  <- marker_map$marker[i]
  mg <- marker_map$marker_group_id[i]

  aa_path <- file.path(OUT_MARKERS, g, m, str_c(mg, "_protein.fasta"))
  nt_path <- file.path(OUT_MARKERS, g, m, str_c(mg, "_nucleotid.fasta"))

  if (!file.exists(aa_path) || !file.exists(nt_path)) {
    next
  }

  AA <- Biostrings::readAAStringSet(aa_path)
  DNA <- Biostrings::readDNAStringSet(nt_path)

  # vsearch
  AA_vs <- AA
  DNA_vs <- DNA

  names(AA_vs) <- sapply(str_split(names(AA_vs), " "), `[`, 1)
  names(DNA_vs) <- sapply(str_split(names(DNA_vs), " "), `[`, 1)

  Biostrings::writeXStringSet(
    AA_vs,
    file.path(outdir, "vsearch", str_c(mg, "_protein.fasta.gz"))
  )

  Biostrings::writeXStringSet(
    DNA_vs,
    file.path(outdir, "vsearch", str_c(mg, "_nucleotid.fasta.gz"))
  )

  # dada2 assignTaxonomy
  AA_dada_gen <- AA
  DNA_dada_gen <- DNA

  taxo_AA <- get_taxo_from_names(AA_dada_gen)
  taxo_DNA <- get_taxo_from_names(DNA_dada_gen)

  names(AA_dada_gen) <- str_c("tax=", taxo_AA)
  names(DNA_dada_gen) <- str_c("tax=", taxo_DNA)

  Biostrings::writeXStringSet(
    AA_dada_gen,
    file.path(outdir, "dada2", str_c(mg, "_train_protein.fasta.gz"))
  )

  Biostrings::writeXStringSet(
    DNA_dada_gen,
    file.path(outdir, "dada2", str_c(mg, "_train_nucleotid.fasta.gz"))
  )

  # dada2 addSpecies
  AA_species <- AA
  DNA_species <- DNA

  split_AA <- str_split(taxo_AA, ";")
  split_DNA <- str_split(taxo_DNA, ";")

  genus_AA <- sapply(split_AA, function(x) ifelse(length(x) >= 6, x[6], "NA"))
  species_AA <- sapply(split_AA, function(x) ifelse(length(x) >= 7, x[7], "NA"))

  genus_DNA <- sapply(split_DNA, function(x) ifelse(length(x) >= 6, x[6], "NA"))
  species_DNA <- sapply(split_DNA, function(x) ifelse(length(x) >= 7, x[7], "NA"))

  names(AA_species) <- str_c(
    sapply(str_split(names(AA_species), " "), `[`, 1),
    genus_AA,
    species_AA,
    sep = " "
  )

  names(DNA_species) <- str_c(
    sapply(str_split(names(DNA_species), " "), `[`, 1),
    genus_DNA,
    species_DNA,
    sep = " "
  )

  Biostrings::writeXStringSet(
    AA_species,
    file.path(outdir, "dada2", str_c(mg, "_species_protein.fasta.gz"))
  )

  Biostrings::writeXStringSet(
    DNA_species,
    file.path(outdir, "dada2", str_c(mg, "_species_nucleotid.fasta.gz"))
  )
}
