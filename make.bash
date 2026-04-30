#!/usr/bin/env bash

# Generate VirMarkDB database 
    # Generate a manifest list of virus
    module load r/4.4.1
    Rscript scripts/analysis/01_generate_manifest.R
    # Search ORF in viral genomes
    bash scripts/analysis/02_prodigal.bash
    
    # Generate HMM profile from references protein
        # Make an alignment of references sequences
        bash scripts/analysis/03_align_markers.bash
        # Build hmm profile
        bash scripts/analysis/04_hmmbuild.bash
        # Search protein with HMM
        bash scripts/analysis/05_hmmsearch.bash
    
    # Generate output files
    Rscript scripts/analysis/06_final_results.R

# Generate a pool of potentential protein from ncbi
    # Load bamford protein from the local nr database of bamfordvirae
    sbatch scripts/analysis/07_scrap_viral_protein_ncbi.bash
    #  Filter out mis annotated protein with an auto alignment
    sbatch scripts/analysis/08_align_pool_vs_pool.bash
    module load r/4.4.1
    Rscript scripts/analysis/09_results_align_pool_vs_pool.R

# Test VirMarkDB-GV with algniment of a pool of viral protein
    # Blast with mmseqs each marker from VirMarkDBgv to each pool prot
    sbatch scripts/analysis/10_align_pool_VirMarkDB.bash
    # Analyze the results
    module load r/4.4.1
    Rscript scripts/analysis/11_analyse_results_align_pool_VirMarkDB.R
    # Clusterise missing protein 
    bash scripts/analysis/12_cluterise_missing_prot.bash
