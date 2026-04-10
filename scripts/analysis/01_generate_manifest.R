library(tidyverse)

# Config
ictv_xlsx <- "input/ICTV_tables/VMR_MSL40.v2.20251013.xlsx"
manual_gen_ncbi <- "input/manual_genomes_ncbi.tsv"
manual_gen_priv <- "input/manual_genomes_private.tsv"
out_tsv <- "output/manifest_genomes.tsv"
genomes_private_path <- "input/private_genomes"
out_genomes <- "output/genomes"


dir.create("output", showWarnings = FALSE)


# Load ICTV data
ictv <- readxl::read_excel(ictv_xlsx, sheet = 2) %>%
  filter(Kingdom == "Bamfordvirae") %>%
  rename(
    virus_id = `Virus GENBANK accession`,
    Virus_names = "Virus name(s)",
    Virus_names_abrv = "Virus name abbreviation(s)",
    Host_source = "Host source"
  ) %>%
  select(
    virus_id, Kingdom, Phylum, Class, Order, Family, Genus, Species,
    Virus_names, Virus_names_abrv, Host_source, ICTV_ID
  ) %>%
  separate_rows(virus_id, sep = "; ") %>%
  mutate(source = "ICTV")


# Load manual genomes on ncbi
if (file.exists(manual_gen_ncbi)) {
  manual <- read_tsv(manual_gen_ncbi, show_col_types = FALSE) %>%
    mutate(source = "MANUAL")

  genomes <- bind_rows(ictv, manual)
} else {
  genomes <- ictv
}

#Check duplication of virus id if exist
genomes %>%
  group_by(virus_id) %>%
  filter(n() > 1) %>%
  pull(virus_id)

# Export
write_tsv(genomes, out_tsv)

# Download genomes by accession
source("scripts/utils/down_ncbi.R") 

# Merge genomes sources
private_files <- list.files(genomes_private_path, full.names = TRUE)


if (length(private_files) > 0) {
  for (file in private_files) {
    file_name <- basename(file) 
    file.copy(from = file, to = file.path(out_genomes, file_name), overwrite = FALSE)
  }
}




# Manage references protein
source("scripts/utils/prepare_ref_protein.R")