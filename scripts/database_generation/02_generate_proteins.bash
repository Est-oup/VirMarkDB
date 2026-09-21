#!/usr/bin/env bash

module load prodigal/2.6.3
module load emboss/6.6.0


MANIFEST="output/config/manifest_genomes.tsv"
GENOME_DIR="output/genomes"
OUT_DIR="output/proteins"

mkdir -p "$OUT_DIR"

awk 'NR>1' "$MANIFEST" | while IFS=$'\t' read -r virus_id Kingdom Phylum Class Order Family Genus Species Virus_names Virus_names_abrv Host_source ICTV_ID Source genome_type
do

    genome="${GENOME_DIR}/${virus_id}.fna"
    protein_out="${OUT_DIR}/${virus_id}.faa"

    if [[ ! -f "$genome" ]]; then
        echo "Missing genome: $genome"
        continue
    fi

    if [[ -s "$protein_out" ]]; then
        echo "Already done: $virus_id"
        continue
    fi

    case "$genome_type" in

        DNA)

            echo "[DNA] Prodigal : $virus_id"

            prodigal \
                -i "$genome" \
                -a "$protein_out" \
                -p meta

            ;;

        RNA)

            echo "[RNA] Transeq : $virus_id"

            transeq \
                -sequence "$genome" \
                -frame 6 \
                -outseq "$protein_out" \
                -auto

            ;;

        *)

            echo "Unknown genome type: $virus_id -> $genome_type"

            ;;

    esac

done