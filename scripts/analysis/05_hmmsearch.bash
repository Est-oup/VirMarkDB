#!/usr/bin/env bash

module load hmmer/3.3.2

FAA_DIR="output/orfs"
HMM_DIR="output/hmm/hmms"
OUT="output/hmm/search"

mkdir -p "$OUT"

for faa in "$FAA_DIR"/*.faa; do
  genome_id="$(basename "$faa" .faa)"

  for hmm in "$HMM_DIR"/*.hmm; do
    marker="$(basename "$hmm" .hmm)"

    hmmsearch -E 0.00001 --tblout "$OUT/${genome_id}__${marker}.tbl" "$hmm" "$faa" > /dev/null
    # hmmsearch --tblout "$OUT/${genome_id}__${marker}.tbl" "$hmm" "$faa" > /dev/null

  done
done
