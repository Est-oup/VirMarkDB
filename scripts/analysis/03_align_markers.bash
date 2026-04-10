#!/usr/bin/env bash

module load mafft/7.525
module load trimal/1.5.0

IN="output/references_protein"
OUT1="output/hmm/aln"
OUT2="output/hmm/aln_filt"

mkdir -p "$OUT1" "$OUT2"

for fasta in "$IN"/*.fasta; do
  stem="$(basename "$fasta" .fasta)"

  mafft --auto "$fasta" > "${OUT1}/${stem}.aln"
  trimal -in "${OUT1}/${stem}.aln" -out "${OUT2}/${stem}.aln" -gt 0.1
done

echo "OK"
echo "MAFFT alignments : $OUT1"
echo "trimAl filtered  : $OUT2"