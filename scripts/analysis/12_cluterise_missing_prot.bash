#!/usr/bin/env bash
#SBATCH --job-name=cdhit_missing_prot
#SBATCH --partition=fast
#SBATCH --cpus-per-task=16
#SBATCH --output=cdhit_missing_prot.out

module load cd-hit/4.8.1

in_dir="output/benchmark/alignment_vmd_pool/missing_prot"
out_dir="output/benchmark/alignment_vmd_pool/missing_prot/clustering"

mkdir -p "${out_dir}"

for fasta in "${in_dir}"/*.fasta; do

    marker=$(basename "${fasta}" .fasta)

    echo "Input : ${fasta}"
    echo "Marker: ${marker}"

    cd-hit \
        -i "${fasta}" \
        -o "${out_dir}/${marker}/${marker}_cdhit95.fasta" \
        -c 0.95 \
        -n 5 \
        -d 0 \
        -M 0 \
        -T "${SLURM_CPUS_PER_TASK}" \
        > "${out_dir}/${marker}/${marker}_cdhit95.log" 2>&1

done