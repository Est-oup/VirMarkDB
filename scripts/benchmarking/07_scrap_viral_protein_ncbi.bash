#!/usr/bin/env bash
#SBATCH --job-name=scrap_viral_prot
#SBATCH --partition=fast
#SBATCH --cpus-per-task=50
#SBATCH --output=scrap_viral_prot.out

module load blast/2.16.0

OUTDIR="output/benchmark/pool_protein/pool_protein_raw"
TMPFILE="output/benchmark/pool_protein/pool_protein_raw/bamford_titles.tsv"

mkdir -p "$OUTDIR"

# SQLite local dowloaded from NCBI Taxonomy FTP on 2026-05-04
LOCAL_TAXDB_DIR="input/ncbi_taxonomy"
# mkdir -p "$LOCAL_TAXDB_DIR"
# wget -O "${LOCAL_TAXDB_DIR}/taxdb.tar.gz" "https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz"
# tar -xzf "${LOCAL_TAXDB_DIR}/taxdb.tar.gz" -C "$LOCAL_TAXDB_DIR"

# NR database shared
SHARED_NR_DIR="/shared/bank/nr/nr_2025-07-20/blast"

# BLAST search on both folder
export BLASTDB="${LOCAL_TAXDB_DIR}:${SHARED_NR_DIR}"

DB="nr"
TAXID="2732005" # taxid of Bamfordvirae

# extract protein name and title of Bamfordvirae
blastdbcmd \
  -db "$DB" \
  -dbtype prot \
  -taxids "$TAXID" \
  -target_only \
  -outfmt "%a---%t" \
  > "$TMPFILE"

# function to filter protein based on the title
fetch_marker() {
  local marker="$1"
  local pattern="$2"
  local exclude='partial|hypothetical|putative|(^|[^A-Za-z])MAG([^A-Za-z]|$)'

  awk -F '---' -v IGNORECASE=1 -v pat="$pattern" -v excl="$exclude" '
    ($2 ~ pat) && ($2 !~ excl) { print $1 }
  ' "$TMPFILE" \
  | sort -u \
  | blastdbcmd \
      -db "$DB" \
      -dbtype prot \
      -entry_batch - \
      -out "${OUTDIR}/${marker}.fasta"
}

# fetch markers
fetch_marker "majorcapsid_nucleocytoviricota" 'major capsid protein|major capsid'
fetch_marker "dnapol_nucleocytoviricota"      'DNA polymerase B|B-family DNA polymerase|family B DNA polymerase|(^|[^A-Za-z])PolB([^A-Za-z]|$)|(^|[^A-Za-z])pPolB([^A-Za-z]|$)'
fetch_marker "atpase_nucleocytoviricota"      'packaging ATPase|DNA packaging ATPase'
fetch_marker "primase_nucleocytoviricota"     '(^|[^A-Za-z])primase([^A-Za-z]|$)|DNA primase'
fetch_marker "vltf3_nucleocytoviricota"       'late transcription factor 3|(^|[^A-Za-z])VLTF3([^A-Za-z]|$)'
fetch_marker "tf2s_nucleocytoviricota"        'transcription elongation factor SII|(^|[^A-Za-z])TFIIS([^A-Za-z]|$)|transcription factor S-II'
fetch_marker "rnapol1_nucleocytoviricota"     'DNA-directed RNA polymerase subunit 1|RNA polymerase subunit 1'
fetch_marker "rnapol2_nucleocytoviricota"     'DNA-directed RNA polymerase subunit 2|RNA polymerase subunit 2'

# final count
echo "capsid   $(grep -c '^>' "$OUTDIR/majorcapsid_nucleocytoviricota.fasta" || true)"
echo "dnapol   $(grep -c '^>' "$OUTDIR/dnapol_nucleocytoviricota.fasta" || true)"
echo "atpase   $(grep -c '^>' "$OUTDIR/atpase_nucleocytoviricota.fasta" || true)"
echo "primase  $(grep -c '^>' "$OUTDIR/primase_nucleocytoviricota.fasta" || true)"
echo "vltf3    $(grep -c '^>' "$OUTDIR/vltf3_nucleocytoviricota.fasta" || true)"
echo "tf2s     $(grep -c '^>' "$OUTDIR/tf2s_nucleocytoviricota.fasta" || true)"
echo "rnapol1  $(grep -c '^>' "$OUTDIR/rnapol1_nucleocytoviricota.fasta" || true)"
echo "rnapol2  $(grep -c '^>' "$OUTDIR/rnapol2_nucleocytoviricota.fasta" || true)"