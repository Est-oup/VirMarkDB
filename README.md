# VMD — Viral Marker Database (for Bamfordvirae dsDNA Viruses)

**VMD** is a reference database of **marker-gene sequences** (**protein + nucleotide**) for dsDNA viruses within **Bamfordvirae**.  
It is intended for **taxonomic assignment**, **marker-based phylogeny**, and more broadly for workflows that require curated reference sets of conserved viral proteins.

This repository documents **how the database was generated** and **how it was refined through a benchmark step** designed to detect missing or underrepresented diversity.  
The final curated release of the database is meant to be distributed separately through a **public archive** (for example Zenodo).

## Public archive

- **Archive / Zenodo:** `[LINK]`
- **DOI:** `[LINK]`

---

## 1) What is inside?

### 1.1 Taxonomic scope

- **Kingdom:** `Bamfordvirae`
- Taxonomy from **Kingdom** to **Species** is derived from the ICTV resources used during the build.

### 1.2 Marker genes included

The database currently includes the following marker sets:

| Marker | Description |
|---|---|
| `majorcapsid` | Major capsid protein (MCP) / capsid morphogenesis |
| `dnapol` | DNA polymerase B / replication |
| `atpase` | Packaging ATPase / virion assembly |
| `primase` | Primase / replication initiation |
| `rnapol1` | RNA polymerase subunit 1 |
| `rnapol2` | RNA polymerase subunit 2 |
| `tf2s` | Transcription factor S-II / TFIIS-like marker |
| `vltf3` | Viral late transcription factor 3 |

### 1.3 Typical use cases

The database can be used for:

- **marker phylogenies** built from curated reference sets,
- **reference FASTA collections** for BLAST, HMM, MMseqs2 or DIAMOND-based workflows,
- **presence/absence summaries** of markers across Bamfordvirae genomes,
- **metadata-aware analyses** relying on ICTV taxonomy and accession-linked provenance.

---

## 2) Repository layout

The workflow currently relies on the following structure:

```bash
.
├── input/
│   ├── ICTV_tables/
│   ├── manual_genomes_ncbi.tsv
│   ├── manual_reference_protein/
│   ├── ncbi.creds
│   └── ncbi_taxonomy/
├── output/
│   ├── benchmark/
│   ├── genomes/
│   ├── hmm/
│   ├── manifest_genomes.tsv
│   ├── orfs/
│   ├── reference_sources/
│   ├── references_protein/
│   └── VMD-database/
│       ├── export_format/
│       ├── markers/
│       └── virus_informations/
├── scripts/
│   ├── analysis/
│   └── utils/
├── make.bash
└── README.md
```

### 2.1 Final database location

The final exported database is written to:

```bash
output/VMD-database/
```

It contains three main components:

- `markers/`
- `virus_informations/`
- `export_format/`

### 2.2 Marker folders

Each marker folder under:

```bash
output/VMD-database/markers/<marker>/
```

contains:

- `<marker>_protein.fasta`
- `<marker>_nucleotid.fasta`
- `<marker>_virus_orf_description.tsv`

### 2.3 Global tables

Global summary tables are written to:

```bash
output/VMD-database/virus_informations/
```

and currently include:

- `virus_compo_taxo.tsv`
- `virus_metadata.tsv`

---

## 3) Keys, headers, and traceability

### 3.1 Primary key: `virus_id`

The central key used throughout the workflow is `virus_id`.

In the current build logic, `virus_id` corresponds to the accession identifier from ncbi and propagated through downstream analyses.  
It is used to join:

- ICTV-derived metadata,
- genome files,
- predicted ORFs,
- HMM search results,
- and final database exports.

### 3.2 FASTA headers

In the current export logic, marker FASTA headers are built from:

- `virus_id`
- followed by the taxonomic lineage from `Kingdom` to `Species`

This makes exported marker sequences directly traceable to their taxonomic context.

---

## 4) Final database content

## 4.1 `virus_compo_taxo.tsv`

This is the main wide-format summary table of the exported database.

It includes:

- `virus_id`
- `Virus_names`
- `ICTV_ID`
- one column per marker
- taxonomy columns from `Kingdom` to `Species`

For each marker column, the stored value corresponds to the **retained ORF name** selected as the best hit for that virus.  
An `NA` means that no hit was retained for that marker in the current build.

### 4.2 `virus_metadata.tsv`

This table contains the metadata associated with viruses retained in the final exported database.

It includes:

- `virus_id`
- `ICTV_ID`
- `Virus_names`
- `Virus_names_abrv`
- `Host_source`
- `Origin_source`
- taxonomy columns from `Kingdom` to `Species`

`Origin_source` indicates how the genome entered the workflow, for example through ICTV-based accession retrieval or manual additions.

### 4.3 Marker-wise folders

For each marker, the database exports:

- a **protein FASTA** of retained hits,
- a **nucleotide FASTA** of retained hits,
- a **TSV table** describing the selected ORFs.

The per-marker description table includes:

- `virus_id`
- `marker`
- `orf_name`
- `evalue`
- `score`
- `prodigal_start`
- `prodigal_end`
- `Species`
- `Virus_names`

These files provide the sequence layer and the minimal hit-level annotation needed for downstream analyses.

---

## 5) Provenance and source data

This build is derived from:

### ICTV resources
- Master Species List: **MSL40.v2** (`ICTV_Master_Species_List_2024_MSL40.v2.xlsx`)  
- Virus Metadata Resource: **VMR aligned to MSL40.v2** (`VMR_MSL40.v2.20251013.xlsx`) 

### Additional inputs
- `input/manual_genomes_ncbi.tsv` for optional manually added NCBI genomes
- `input/manual_reference_protein/` for manually curated protein references
- local NCBI-related resources in `input/ncbi_taxonomy/`

### Reference protein sources

The HMM profiles are built from marker-specific reference protein sets stored in:

```bash
output/references_protein/
```

These reference sets are part of the documented generation process and can be complemented by curated manual additions.

---

## 6) How the database was built

The build process can be summarised in two major phases:

1. **initial database construction**
2. **benchmark-guided refinement**

This distinction is important: the final database is not only the result of a one-pass marker detection workflow, but also of an additional evaluation phase designed to reveal what the initial database was still missing.

### 6.1 Initial construction

The first phase consisted in generating a structured marker database directly from Bamfordvirae genomes.

Main steps:

1. **Manifest generation**  
   Bamfordvirae entries are extracted from ICTV tables, accession information is normalised, and a consolidated manifest is written to:

   ```bash
   output/manifest_genomes.tsv
   ```

2. **Genome retrieval**  
   Genome sequences associated with the manifest are downloaded and stored in:

   ```bash
   output/genomes/
   ```

3. **ORF prediction**  
   ORFs are predicted from each genome using **Prodigal** in metagenomic mode (`-p meta`), generating:

   - protein FASTA files,
   - nucleotide FASTA files,
   - GFF annotation files.

4. **Reference alignment and HMM construction**  
   Marker reference proteins are aligned with **MAFFT**, filtered with **trimAl**, and converted into marker-specific HMM profiles with **HMMER hmmbuild**.

5. **Marker detection in predicted ORFs**  
   Predicted ORF proteomes are searched with **hmmsearch** using the generated marker HMM profiles.

6. **Best-hit selection and final export**  
   For each `(virus_id, marker)` pair, the best hit is selected using HMM search ranking and exported as part of the final database structure.

### 6.2 Benchmark-guided refinement

A second phase was used to evaluate the coverage of the initial database against a broader Bamfordvirae protein pool extracted from a local **NR** database.

This benchmark phase had a practical goal:  
**detect the holes in the database** in other words, identify which parts of marker diversity were still missing or insufficiently represented in the initial VMD-GV build.

Main steps:

1. **Extraction of a broader Bamfordvirae protein pool from NR**  
   A local BLAST-enabled NR database is queried using Bamfordvirae taxonomic restriction, and marker candidates are extracted through title-based filtering.

2. **Self-comparison of each marker pool**  
   Each marker pool is aligned against itself with **MMseqs2** to identify isolated, weakly connected, or suspicious sequences.

3. **Filtering of the benchmark pool**  
   Proteins that do not show satisfactory similarity to other members of the same marker pool are removed, producing a cleaner benchmark set.

4. **Comparison of the filtered pool against the current VMD-GV database**  
   The filtered external pool is then aligned against the marker sequences already present in the database.

5. **Identification of missing proteins**  
   Proteins that fail to match the current database are extracted as potentially missing diversity.

6. **Clustering of missing proteins**  
   These missing proteins are clustered with **CD-HIT** to reduce redundancy and facilitate downstream review.

### 6.3 Why the benchmark matters

This second phase is central to the history of the database.

It was used to:

- evaluate how much external Bamfordvirae marker diversity was already covered,
- detect proteins not captured by the initial database,
- identify underrepresented regions of marker space,
- and improve the reference space used to generate the final VMD-GV release.

In practice, this benchmark step served to **improve the HMM profiles indirectly** by revealing what was missing and therefore what needed to be better represented.

---

## 7) Workflow summary

The overall execution logic is summarised in:

```bash
make.bash
```

This script records the sequence of operations used during database generation and benchmark evaluation:

Database generation:
1. manifest generation,
2. ORF prediction,
3. marker reference alignment,
4. HMM construction,
5. HMM search on predicted ORFs,
6. export of the initial VMD-GV database,

Database benchmarking:
7. extraction of a Bamfordvirae protein pool from NR,
8. self-alignment and filtering of that pool,
9. comparison of the filtered pool against the database,
10. analysis of missing proteins,
11. clustering of missing candidates.

---

## 8) Software used

The workflow relies on the following tools:

- **R**
- **Prodigal**
- **MAFFT**
- **trimAl**
- **HMMER**
- **BLAST+**
- **MMseqs2**
- **CD-HIT**

Examples of module versions visible in the workflow include:

- `r/4.4.1`
- `prodigal/2.6.3`
- `mafft/7.525`
- `trimal/1.5.0`
- `hmmer/3.3.2`
- `blast/2.16.0`
- `mmseqs2/15.6f452`
- `cd-hit/4.8.1`

---

## 9) Notes and limitations

- **ORF prediction** was performed with Prodigal, which is practical and robust, but viral genomes can still present difficult coding configurations.
- **HMM profile performance** depends strongly on the diversity and quality of the seed reference sets.
- The exported database keeps **one retained hit per virus and per marker** in the final main exports.
- `NA` values in marker presence tables should be interpreted as **not detected under the current workflow and thresholds**, not automatically as true biological absence.
- The benchmark filtering strategy is **empirical** and was designed as a practical way to clean large external protein pools before coverage assessment.

---

## 10) Recommended citation and license

To be completed when the public archive is deposited.

Suggested placeholders:

- **Database citation:** `[CITATION]`
- **License:** `[LICENSE]`

---

## 11) References

### ICTV resources
- ICTV Master Species Lists (MSL)
- ICTV Virus Metadata Resource (VMR)

### Core tools used in the workflow
- Prodigal
- MAFFT
- trimAl
- HMMER
- BLAST+
- MMseqs2
- CD-HIT

### Reference datasets and literature
- Guglielmini et al. (2019), on large and giant eukaryotic dsDNA viruses
- Associated Zenodo supplementary data deposit used as one of the reference sources (Zenodo)[https://zenodo.org/records/3368642/files/Additional%20data.zip?download=1]

---

## Contact / issues

This repository is intended to document the generation and refinement of the VMD-GV database.

Potential uses of the issue tracker include:

- suspicious marker assignments,
- missing taxa or poorly represented lineages,
- requests for additional marker sets,
- questions about workflow provenance or reproducibility.

Feel free to suggest something it will be very helpful !