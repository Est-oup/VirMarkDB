library(tidyverse)
library(fs)

outdir <- "output/VMD-database/export_format/"

# Build taxonomy file for vearch

taxonomy_all <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

taxo_vsearch <- virus_compo_taxo %>%
  select(virus_id, any_of(taxonomy_all)) %>%
  mutate(reference_id = virus_id) %>%
  unite(taxonomy, any_of(taxonomy_all), sep = ";") %>%
  select(reference_id, taxonomy)

# export taxonomy
dir_create(str_c(outdir,"/vsearch"))
write_tsv(taxo_vsearch, str_c(outdir,"/vsearch/taxonomy.tsv"))

# export sequences for dada2 and vsearch 
dir_create(str_c(outdir,"/dada2"))

for (m in markers){
    AA <- Biostrings::readAAStringSet(str_c("output/VMD-database/markers/",m,"/",m,"_protein.fasta"))
    DNA <- Biostrings::readDNAStringSet(str_c("output/VMD-database/markers/",m,"/",m,"_nucleotid.fasta"))

    #vsearch
    AA_vs <- AA
    DNA_vs <- DNA

    names(AA_vs) <- sapply(str_split(names(AA_vs)," "), `[`, 1)
    names(DNA_vs) <- sapply(str_split(names(DNA_vs)," "), `[`, 1)

    Biostrings::writeXStringSet(AA_vs, str_c(outdir,"/vsearch/",m,"_protein.fasta.gz"))
    Biostrings::writeXStringSet(DNA_vs, str_c(outdir,"/vsearch/",m,"_nucleotid.fasta.gz"))
    
    # dada2 : assignTaxonomy
    AA_dada_gen  <- AA
    DNA_dada_gen <- DNA

    taxo_AA  <- sapply(str_split(names(AA_dada_gen),  " "), `[`, 2)   
    taxo_DNA <- sapply(str_split(names(DNA_dada_gen), " "), `[`, 2)

    names(AA_dada_gen)  <- str_c("tax=", taxo_AA)
    names(DNA_dada_gen) <- str_c("tax=", taxo_DNA)

    # dada2 : toSpecies
    AA_species  <- AA
    DNA_species <- DNA

    split_AA  <- str_split(taxo_AA,  ";")
    split_DNA <- str_split(taxo_DNA, ";")

    genus_AA   <- sapply(split_AA,  `[`, 6)
    species_AA <- sapply(split_AA,  `[`, 7)

    genus_DNA   <- sapply(split_DNA, `[`, 6)
    species_DNA <- sapply(split_DNA, `[`, 7)

    names(AA_species)  <- str_c(sapply(str_split(names(AA_species)," "), `[`, 1),genus_AA,  species_AA,  sep = " ")
    names(DNA_species) <- str_c(sapply(str_split(names(DNA_species)," "), `[`, 1),genus_DNA, species_DNA, sep = " ")

    Biostrings::writeXStringSet(AA_dada_gen,  str_c(outdir, "/dada2/", m, "_train_protein.fasta.gz"))
    Biostrings::writeXStringSet(DNA_dada_gen, str_c(outdir, "/dada2/", m, "_train_nucleotid.fasta.gz"))
    Biostrings::writeXStringSet(AA_species,  str_c(outdir, "/dada2/", m, "_species_protein.fasta.gz"))
    Biostrings::writeXStringSet(DNA_species, str_c(outdir, "/dada2/", m, "_species_nucleotid.fasta.gz"))
}


