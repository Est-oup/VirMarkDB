#!/usr/bin/env bash
#SBATCH --job-name=align_vmd_pool
#SBATCH --partition=fast
#SBATCH --cpus-per-task=50
#SBATCH --output=align_vmd_pool.out

set -euo pipefail

module load mmseqs2/15.6f452

mkdir -p output/benchmark/alignment_vmd_pool/alignment/tmp

vmd_db="output/VMD-database/markers"
pool_dir="output/benchmark/pool_protein/pool_protein_filt"
out_dir="output/benchmark/alignment_vmd_pool/alignment"

for ref in "$vmd_db"/*/*/*_protein.fasta; do
    full_name=$(basename "$ref" _protein.fasta)
    name="${full_name##*_}"
    marker="${full_name%_*}"
    query="${pool_dir}/${marker}.fasta"

    echo "Running marker_group_id: $name"
    echo "Using pool marker: $marker"

    mmseqs easy-search \
        "$query" \
        "$ref" \
        "${out_dir}/${marker}_${name}.m8" \
        output/benchmark/alignment_vmd_pool/alignment/tmp \
        --threads 50
done

