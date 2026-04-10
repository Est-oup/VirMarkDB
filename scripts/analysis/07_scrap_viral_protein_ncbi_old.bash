#!/usr/bin/env bash

module load blast/2.16.0

OUTDIR=output/benchmark/pool_protein/pool_protein_raw
mkdir -p "$OUTDIR"

TAXON='txid2732005[Organism:exp]'
EXCLUDE='NOT (putative[Title] OR partial[Title] OR MAG[Title] OR hypothetical[Title])'

esearch -db protein -query "$TAXON AND (\"major capsid protein\"[Title] OR \"major capsid\"[Title]) $EXCLUDE" \
| efetch -format fasta > "$OUTDIR/majorcapsid.fasta"

esearch -db protein -query "$TAXON AND (\"DNA polymerase B\"[Title] OR \"B-family DNA polymerase\"[Title] OR \"family B DNA polymerase\"[Title] OR PolB[Title] OR pPolB[Title]) $EXCLUDE" \
| efetch -format fasta > "$OUTDIR/dnapol.fasta"

esearch -db protein -query "$TAXON AND (\"packaging ATPase\"[Title] OR \"DNA packaging ATPase\"[Title]) $EXCLUDE" \
| efetch -format fasta > "$OUTDIR/atpase.fasta"

esearch -db protein -query "$TAXON AND (primase[Title] OR \"DNA primase\"[Title]) $EXCLUDE" \
| efetch -format fasta > "$OUTDIR/primase.fasta"

esearch -db protein -query "$TAXON AND (\"late transcription factor 3\"[Title] OR VLTF3[Title]) $EXCLUDE" \
| efetch -format fasta > "$OUTDIR/vltf3.fasta"

esearch -db protein -query "$TAXON AND (\"transcription elongation factor SII\"[Title] OR TFIIS[Title] OR \"transcription factor S-II\"[Title]) $EXCLUDE" \
| efetch -format fasta > "$OUTDIR/tf2s.fasta"

esearch -db protein -query "$TAXON AND (\"DNA-directed RNA polymerase subunit 1\"[Title] OR \"RNA polymerase subunit 1\"[Title]) $EXCLUDE" \
| efetch -format fasta > "$OUTDIR/rnapol1.fasta"

esearch -db protein -query "$TAXON AND (\"DNA-directed RNA polymerase subunit 2\"[Title] OR \"RNA polymerase subunit 2\"[Title]) $EXCLUDE" \
| efetch -format fasta > "$OUTDIR/rnapol2.fasta"

echo "capsid   $(grep -c '^>' "$OUTDIR/majorcapsid.fasta" || true)"
echo "dnapol   $(grep -c '^>' "$OUTDIR/dnapol.fasta" || true)"
echo "atpase   $(grep -c '^>' "$OUTDIR/atpase.fasta" || true)"
echo "primase  $(grep -c '^>' "$OUTDIR/primase.fasta" || true)"
echo "vltf3    $(grep -c '^>' "$OUTDIR/vltf3.fasta" || true)"
echo "tf2s     $(grep -c '^>' "$OUTDIR/tf2s.fasta" || true)"
echo "rnapol1  $(grep -c '^>' "$OUTDIR/rnapol1.fasta" || true)"
echo "rnapol2  $(grep -c '^>' "$OUTDIR/rnapol2.fasta" || true)"