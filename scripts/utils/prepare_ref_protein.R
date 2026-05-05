library(tidyverse)
library(rentrez)

zenodo_dir <- "output/reference_sources/zenodo_raw"
zenodo_files_path <- "output/reference_sources/zenodo_raw/Alignments/Sequences"
manual_ids_protein_path <- "input/manual_reference_protein"
merged_sequences_path <- "output/references_protein"

dir.create(zenodo_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(merged_sequences_path, recursive = TRUE, showWarnings = FALSE)

# Download and unzip Zenodo references
zip_url  <- "https://zenodo.org/records/3368642/files/Additional%20data.zip?download=1"
zip_dest <- file.path(zenodo_dir, "Additional_data.zip")

if (!file.exists(zip_dest)) {
  download.file(zip_url, destfile = zip_dest, mode = "wb")
  unzip(zip_dest, exdir = zenodo_dir)
}


# Load Zenodo protein references
majorcapsid_nucleocytoviricota_zenodo <- Biostrings::readAAStringSet(file.path(zenodo_files_path, "capsid.prt"))
atpase_nucleocytoviricota_zenodo       <- Biostrings::readAAStringSet(file.path(zenodo_files_path, "pATPase.prt"))
dnapol_nucleocytoviricota_zenodo       <- Biostrings::readAAStringSet(file.path(zenodo_files_path, "dnapol.prt"))
primase_nucleocytoviricota_zenodo      <- Biostrings::readAAStringSet(file.path(zenodo_files_path, "primase.prt"))
rnapol1_nucleocytoviricota_zenodo      <- Biostrings::readAAStringSet(file.path(zenodo_files_path, "rnapol1.prt"))
rnapol2_nucleocytoviricota_zenodo      <- Biostrings::readAAStringSet(file.path(zenodo_files_path, "rnapol2.prt"))
tf2s_nucleocytoviricota_zenodo         <- Biostrings::readAAStringSet(file.path(zenodo_files_path, "tf2s.prt"))
vltf3_nucleocytoviricota_zenodo        <- Biostrings::readAAStringSet(file.path(zenodo_files_path, "vltf3.prt"))

atpase_polinto_zenodo                  <- Biostrings::readAAStringSet(file.path(zenodo_files_path, "pATPase-polinto.prt"))
majorcapsid_polinto_zenodo             <- Biostrings::readAAStringSet(file.path(zenodo_files_path, "capsid-polinto.prt"))

# Load manual accessions and fetch NCBI protein sequences
manual_markers <- list.files(
  manual_ids_protein_path,
  pattern = "\\.tsv$",
  full.names = FALSE
) %>%
  tools::file_path_sans_ext()




for (marker in manual_markers) {
  print(marker)
  manual_file <- str_c(manual_ids_protein_path, "/", marker, ".tsv")

  if (file.info(manual_file)$size == 0) {
    assign(str_c(marker, "_manual"), tibble(
      accession = character(),
      description = character()
    ))
    assign(str_c(marker, "_seq"), Biostrings::AAStringSet())
    next
  }
  
  x <- manual_file |>
    read_tsv(col_names = c("accession", "description"),show_col_types = FALSE)

  assign(str_c(marker, "_manual"), x)
  
  if (nrow(x) == 0) {
    assign(str_c(marker, "_seq"), Biostrings::AAStringSet())
    next
  }
  
  fasta_txt <- entrez_fetch(
    db = "protein",
    id = x$accession,
    rettype = "fasta",
    retmode = "text"
  )
  
  tmp_fasta <- tempfile(fileext = ".fasta")
  writeLines(fasta_txt, tmp_fasta)
  
  seq_obj <- Biostrings::readAAStringSet(tmp_fasta)
  
  assign(str_c(marker, "_seq"), seq_obj)
}

# Merge Zenodo + manual sequences, remove duplicates, write FASTA

all_markers <- c("majorcapsid_nucleocytoviricota",
"atpase_nucleocytoviricota",      
"dnapol_nucleocytoviricota",
"primase_nucleocytoviricota",
"rnapol1_nucleocytoviricota",
"rnapol2_nucleocytoviricota",
"tf2s_nucleocytoviricota",
"vltf3_nucleocytoviricota") 

for (marker in all_markers) {
  
  zenodo_name <- str_c(marker, "_zenodo")
  manual_name <- str_c(marker, "_seq")
  
  zenodo_seq <- if (exists(zenodo_name, inherits = FALSE)) get(zenodo_name) else Biostrings::AAStringSet()
  manual_seq <- if (exists(manual_name, inherits = FALSE)) get(manual_name) else Biostrings::AAStringSet()
  
  merged_seq <- c(zenodo_seq, manual_seq)
  merged_seq <- merged_seq[!duplicated(as.character(merged_seq))]
  
  assign(marker, merged_seq)
  
  Biostrings::writeXStringSet(
    merged_seq,
    filepath = file.path(merged_sequences_path, str_c(marker, ".fasta")),
    format = "fasta"
  )
  
  cat(marker, ":", length(merged_seq), "sequences written\n")
}
