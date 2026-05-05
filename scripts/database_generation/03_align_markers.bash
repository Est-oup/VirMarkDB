#!/usr/bin/env bash

module load mafft/7.525
module load trimal/1.5.0

IN="output/references_protein"
OUT1="output/hmm/aln"
OUT2="output/hmm/aln_filt"

mkdir -p "$OUT1" "$OUT2"

for fasta in "$IN"/*.fasta; do
  stem="$(basename "$fasta" .fasta)"

  mafft --ep 0 --genafpair --maxiterate 1000 "$fasta" > "${OUT1}/${stem}.aln"

done

echo "OK"
echo "MAFFT alignments : $OUT1"
echo "trimAl filtered  : $OUT2"