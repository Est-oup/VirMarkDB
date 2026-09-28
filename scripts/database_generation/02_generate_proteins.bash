#!/usr/bin/env bash

module load prodigal/2.6.3
module load emboss/6.6.0

GENOME_DIR="output/genomes"
OUT_DIR="output/proteins"
MANIFEST="output/config/manifest_genomes.tsv"

mkdir -p "$OUT_DIR"

for genome in "$GENOME_DIR"/*.fna; do
    id=$(basename "$genome" .fna)
    genome_type=$(awk -v id="$id" '$1 == id {print $NF; exit}' "$MANIFEST")

    if [[ "$genome_type" == "DNA" ]]; then
        faa="$OUT_DIR/${id}.faa"
        fna="$OUT_DIR/${id}.fna"
        gff="$OUT_DIR/${id}.gff"

        if [[ -s "$faa" && -s "$fna" && -s "$gff" ]]; then
            echo "[SKIP] Already done: $id"
            continue
        fi

        echo "[DNA] Prodigal: $id"
        prodigal -i "$genome" -a "$faa" -d "$fna" -f gff -o "$gff" -p meta

    elif [[ "$genome_type" == "RNA" ]]; then
        faa="$OUT_DIR/${id}.faa"

        if [[ -s "$faa" ]]; then
            echo "[SKIP] Already done: $id"
            continue
        fi

        echo "[RNA] Transeq: $id"
        transeq -sequence "$genome" -frame 6 -outseq "$faa" -auto

    else
        echo "[ERROR] Unknown genome type for $id: $genome_type" >&2
    fi
done