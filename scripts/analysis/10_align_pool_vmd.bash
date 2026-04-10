#!/usr/bin/env bash
#SBATCH --job-name=align_vmd_pool
#SBATCH --partition=fast
#SBATCH --cpus-per-task=50
#SBATCH --output=align_vmd_pool.out

module load mmseqs2/15.6f452

mkdir -p output/benchmark/alignment_vmd_pool/alignment/tmp

vmd_db="output/VMD-database"

for marker in "$vmd_db"/markers/*; do
    
    name=$(basename "$marker")
    
    echo "Running marker: $name"
    
    mmseqs easy-search \
        "output/benchmark/pool_protein/pool_protein_filt/${name}.fasta" \
        "$marker/${name}_protein.fasta" \
        "output/benchmark/alignment_vmd_pool/alignment/${name}.m8" \
        output/benchmark/alignment_vmd_pool/alignment/tmp \
        --threads 50

done
