#!/usr/bin/env bash

module load blast/2.16.0


# =========================
# A ADAPTER SI BESOIN
# =========================

# Dossier perso où tu veux stocker la taxo NCBI
LOCAL_TAXDB_DIR="input/ncbi_nr"

# Dossier partagé contenant la base nr
# Remplace si besoin par le vrai chemin de ton cluster
SHARED_NR_DIR="/shared/bank/nr/nr_2025-07-20/blast"

# Nom de la base BLAST
DB_NAME="nr"

# TaxID Bamfordvirae
TAXID="2732005"

# =========================
# PREPARATION
# =========================

mkdir -p "${LOCAL_TAXDB_DIR}"
cd "${LOCAL_TAXDB_DIR}"

# Téléchargement de taxdb si absent
if [[ ! -f "taxonomy4blast.sqlite3" ]]; then
    echo "[INFO] Téléchargement de taxdb.tar.gz dans ${LOCAL_TAXDB_DIR}"
    wget -O taxdb.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
    tar -xzf taxdb.tar.gz
fi

# Vérification
if [[ ! -f "taxonomy4blast.sqlite3" ]]; then
    echo "[ERROR] taxonomy4blast.sqlite3 introuvable après extraction"
    exit 1
fi

# On définit BLASTDB avec :
# 1) le dossier perso contenant taxonomy4blast.sqlite3
# 2) le dossier partagé contenant nr
export BLASTDB="${PWD}:${SHARED_NR_DIR}"


blastdbcmd -info -db "${DB_NAME}" -dbtype prot

blastdbcmd \
  -db "${DB_NAME}" \
  -dbtype prot \
  -taxids "${TAXID}" \
  -target_only \
  -outfmt "%a %T %t" | head -n 20

