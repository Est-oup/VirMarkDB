library(tidyverse)
library(rentrez)

# Config
manifest <- "output/config/manifest_genomes.tsv"
outdir   <- "output/genomes"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# NCBI credentials
creds <- read_lines("input/ncbi.creds")
options(entrez_email = creds[1])
options(entrez_key   = creds[2])

# Load manifest and filter out private genomes
genomes <- read_tsv(manifest, show_col_types = FALSE) %>%
  filter(Source != "PRIVATE")

# Function to download genomes
download_genome <- function(virus_id) {
  outfile <- file.path(outdir, paste0(virus_id, ".fna"))

  if (file.exists(outfile)) {
    return()
  }

  seq <- rentrez::entrez_fetch(
    db = "nuccore",
    id = gsub("_partial","",virus_id),
    rettype = "fasta",
    retmode = "text"
  )
  write_lines(seq, outfile)
}

# Apply funtion to download the genomes
genomes %>%
  select(virus_id) %>%
  pwalk(download_genome)
