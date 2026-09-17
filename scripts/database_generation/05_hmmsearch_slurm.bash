#!/usr/bin/env bash
#SBATCH --cpus-per-task=70
#SBATCH --mem=100G
#SBATCH --output=hmm_search_%j.out

module load hmmer/3.3.2

ORF_DIR="output/orfs"
HMM_DIR="output/hmm/hmms"

MAP_MARKERS="input/mapfile.tsv"
MAP_GENOMES="output/config/genome_marker_map.tsv"

OUT_DIR="output/hmm/search"
POOL_DIR="output/hmm/pools"

mkdir -p "$OUT_DIR" "$POOL_DIR"

# Pour chaque marker_group_id présent dans input/mapfile.tsv
tail -n +2 "$MAP_MARKERS" | cut -f5 | sort -u | while read -r marker; do

  echo "[MARKER] $marker"

  hmm="$HMM_DIR/${marker}.hmm"
  pool="$POOL_DIR/${marker}.faa"
  tbl="$OUT_DIR/${marker}.tbl"
  log="$OUT_DIR/${marker}.log"

  # Crée le pool protéique des génomes concernés par ce marker
  > "$pool"

  awk -F'\t' -v marker="$marker" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "virus_id" || $i == "genome_id") genome_col = i
        if ($i == "marker_group_id") marker_col = i
      }
      next
    }

    $marker_col == marker {
      print $genome_col
    }
  ' "$MAP_GENOMES" | sort -u | while read -r genome; do

    faa="$ORF_DIR/${genome}.faa"

    if [[ -s "$faa" ]]; then
      cat "$faa" >> "$pool"
    else
      echo "[WARN] missing ORF file: $faa"
    fi

  done

  echo "[SEARCH] hmmsearch $marker"

  hmmsearch \
    --cpu "${SLURM_CPUS_PER_TASK:-1}" \
    --noali \
    -E 1e-5 \
    --domtblout "$tbl" \
    "$hmm" \
    "$pool" \
    > "$log"

  echo "[OK] $tbl"

done
