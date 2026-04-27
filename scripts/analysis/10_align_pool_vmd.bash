#!/usr/bin/env bash
#SBATCH --job-name=align_vmd_pool
#SBATCH --partition=fast
#SBATCH --cpus-per-task=50
#SBATCH --output=align_vmd_pool.out

module load mmseqs2/15.6f452

mkdir -p output/benchmark/alignment_vmd_pool/alignment/tmp

vmd_db="output/VMD-database/markers"
pool_dir="output/benchmark/pool_protein/pool_protein_filt"
out_dir="output/benchmark/alignment_vmd_pool/alignment"

for ref in "$vmd_db"/*/*/*_protein.fasta; do
    full_name=$(basename "$ref" _protein.fasta)
    query="${pool_dir}/${full_name}.fasta"

    echo "Running marker_group_id: $full_name"
    echo "Using pool marker: $full_name"
    echo "$full_name"
    echo "$ref"

    mmseqs easy-search \
        "$query" \
        "$ref" \
        "${out_dir}/${full_name}.m8" \
        output/benchmark/alignment_vmd_pool/alignment/tmp \
        --threads 50
done

