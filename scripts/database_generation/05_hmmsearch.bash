#!/usr/bin/env bash

module load hmmer/3.3.2

FAA_DIR="output/orfs"
HMM_DIR="output/hmm/hmms"
MAP="output/config/genome_marker_map.tsv"
OUT="output/hmm/search"

mkdir -p "$OUT"

marker_col=$(head -n 1 "$MAP" | tr '\t' '\n' | grep -nx "marker_group_id" | cut -d: -f1)

for faa in "$FAA_DIR"/*.faa; do
  genome_id="$(basename "$faa" .faa)"

  grep "^${genome_id}" "$MAP" | cut -f"$marker_col" | grep -v "^NA$" | sort -u | while read -r marker_group_id; do

    hmm="${HMM_DIR}/${marker_group_id}.hmm"

    echo "Genome: $genome_id | Marker: $marker_group_id"

    hmmsearch -E 0.00001 --domtblout "$OUT/${genome_id}__${marker_group_id}.tbl" "$hmm" "$faa" > /dev/null
  done
done
