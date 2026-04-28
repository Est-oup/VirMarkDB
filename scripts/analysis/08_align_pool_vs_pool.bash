#!/usr/bin/env bash
#SBATCH --job-name=align_vmd_vmd
#SBATCH --partition=fast
#SBATCH --cpus-per-task=50
#SBATCH --output=align_vmd_vmd.out

module load mmseqs2/15.6f452

mkdir -p output/benchmark/pool_protein/auto_alignment/tmp/tmp_aln
mkdir -p output/benchmark/pool_protein/auto_alignment/blast

pool="output/benchmark/pool_protein/pool_protein_raw/"

for marker in "$pool"*.fasta; do
    name=$(basename "$marker" .fasta)

    echo "Running: $name"

    mmseqs easy-search \
        "$marker" \
        "$marker" \
        "output/benchmark/pool_protein/auto_alignment/blast/${name}.m8" \
        output/benchmark/pool_protein/auto_alignment/tmp/tmp_aln \
        --threads 50 
done