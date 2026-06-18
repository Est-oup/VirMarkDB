#!/usr/bin/env Rscript

# Build analysis tables for the VirMarkDB website.

args <- commandArgs(trailingOnly = FALSE)
script_arg <- "--file="
script_path <- sub(script_arg, "", args[grep(script_arg, args)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()
docs_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
repo_root <- normalizePath(file.path(docs_dir, ".."), mustWork = FALSE)
out_dir <- file.path(repo_root, "output", "VirMarkDB")
data_dir <- file.path(docs_dir, "data")

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

write_empty <- function(path, columns) {
  empty <- as.data.frame(setNames(replicate(length(columns), character(0), simplify = FALSE), columns))
  write.table(empty, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

if (!dir.exists(out_dir)) {
  warning("output/VirMarkDB/ was not found. Only template tables were kept.")
  write_empty(file.path(data_dir, "db_file_inventory.tsv"), c("path", "type", "size_bytes"))
  write_empty(file.path(data_dir, "db_marker_group_summary.tsv"), c("marker_group_id", "marker", "group_id", "n_orfs", "n_genomes", "n_families"))
  write_empty(file.path(data_dir, "db_taxonomy_summary.tsv"), c("rank", "n_values"))
  quit(save = "no", status = 0)
}

# File inventory.
files <- list.files(out_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
file_info <- file.info(files)
file_inventory <- data.frame(
  path = sub(paste0(normalizePath(out_dir), "/?"), "", normalizePath(files)),
  type = tools::file_ext(files),
  size_bytes = file_info$size,
  stringsAsFactors = FALSE
)
write.table(file_inventory, file.path(data_dir, "db_file_inventory.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Marker-group ORF summaries.
orf_files <- list.files(out_dir, pattern = "_virus_orf_description\\.tsv$", recursive = TRUE, full.names = TRUE)
marker_summaries <- list()

for (file in orf_files) {
  tab <- tryCatch(read.delim(file, check.names = FALSE), error = function(e) NULL)
  if (is.null(tab) || nrow(tab) == 0) next

  marker <- if ("marker" %in% names(tab)) unique(tab$marker)[1] else NA_character_
  group_id <- if ("group_id" %in% names(tab)) unique(tab$group_id)[1] else NA_character_
  marker_group_id <- sub("_virus_orf_description\\.tsv$", "", basename(file))
  n_genomes <- if ("virus_id" %in% names(tab)) length(unique(tab$virus_id)) else NA_integer_
  n_families <- if ("Family" %in% names(tab)) length(unique(tab$Family)) else NA_integer_

  marker_summaries[[length(marker_summaries) + 1]] <- data.frame(
    marker_group_id = marker_group_id,
    marker = marker,
    group_id = group_id,
    n_orfs = nrow(tab),
    n_genomes = n_genomes,
    n_families = n_families,
    stringsAsFactors = FALSE
  )
}

if (length(marker_summaries) > 0) {
  marker_summary <- do.call(rbind, marker_summaries)
} else {
  marker_summary <- data.frame(marker_group_id = character(), marker = character(), group_id = character(), n_orfs = integer(), n_genomes = integer(), n_families = integer())
}
write.table(marker_summary, file.path(data_dir, "db_marker_group_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Taxonomy summary from metadata.
metadata_file <- file.path(out_dir, "virus_informations", "virus_metadata.tsv")
if (file.exists(metadata_file)) {
  meta <- read.delim(metadata_file, check.names = FALSE)
  ranks <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), names(meta))
  tax_summary <- data.frame(
    rank = ranks,
    n_values = sapply(ranks, function(x) length(unique(meta[[x]][!is.na(meta[[x]]) & meta[[x]] != "NA"]))),
    stringsAsFactors = FALSE
  )
} else {
  tax_summary <- data.frame(rank = character(), n_values = integer())
}
write.table(tax_summary, file.path(data_dir, "db_taxonomy_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Coverage table from virus_compo_taxo.tsv when possible.
compo_file <- file.path(out_dir, "virus_informations", "virus_compo_taxo.tsv")
if (file.exists(compo_file)) {
  compo <- read.delim(compo_file, check.names = FALSE)
  marker_cols <- grep("_orf_names$", names(compo), value = TRUE)
  if (length(marker_cols) > 0 && "Family" %in% names(compo)) {
    markers <- sub("_.*", "", marker_cols)
    families <- sort(unique(compo$Family))
    coverage <- data.frame(Family = families, stringsAsFactors = FALSE)
    for (i in seq_along(marker_cols)) {
      col <- marker_cols[i]
      marker <- markers[i]
      coverage[[marker]] <- sapply(families, function(fam) {
        x <- compo[compo$Family == fam, col]
        sum(!is.na(x) & x != "NA" & x != "")
      })
    }
    totals <- data.frame(Family = "Total genomes", as.list(colSums(coverage[, -1, drop = FALSE])), check.names = FALSE)
    coverage <- rbind(coverage, totals)
    write.table(coverage, file.path(data_dir, "marker_coverage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
