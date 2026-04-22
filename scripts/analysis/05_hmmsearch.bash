#!/usr/bin/env bash

module load hmmer/3.3.2

FAA_DIR="output/orfs"
HMM_DIR="output/hmm/hmms"
MAP="output/config/genome_marker_map.tsv"
OUT="output/hmm/search"

mkdir -p "$OUT"

for faa in "$FAA_DIR"/*.faa; do
  genome_id="$(basename "$faa" .faa)"

  grep "^${genome_id}" "$MAP" | cut -f11 | sort -u | while read -r marker_group_id; do

    if [ -z "$marker_group_id" ]; then
      continue
    fi

    hmm="${HMM_DIR}/${marker_group_id}.hmm"

    if [ ! -f "$hmm" ]; then
      echo "Missing HMM: $hmm"
      continue
    fi

    echo "Genome: $genome_id | Marker: $marker_group_id"

    hmmsearch -E 0.00001 --tblout "$OUT/${genome_id}__${marker_group_id}.tbl" "$hmm" "$faa" > /dev/null
  done
done

