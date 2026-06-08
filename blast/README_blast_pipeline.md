# 🔁 Reiterative BLASTp Pipeline

A two-step R pipeline for performing **reiterative local BLASTp searches** across one or multiple species proteomes, and consolidating the unique hits into a single annotated FASTA file ready for downstream analyses (multiple sequence alignment, phylogenetics, etc.).

---

## Pipeline overview

```
Query protein (FASTA)
        │
        ▼
┌─────────────────────────────┐
│  reiterative_BLASTp_local.R │  ── Step 1: iterative BLASTp search
└─────────────────────────────┘
        │
        ▼
  Per-species BLAST results (.xlsx)
        │
        ▼
┌──────────────────────────────────────┐
│  extract_unique_blast_sequences.R    │  ── Step 2: deduplication & FASTA merge
└──────────────────────────────────────┘
        │
        ▼
  todas_las_secuencias.faa  (all unique hits, species-annotated)
```

---

## Step 1 — `reiterative_BLASTp_local.R`

Runs a local BLASTp between a query protein and a subject proteome, then uses each filtered hit as a new query to expand the search iteratively — similar to PSI-BLAST logic but fully controlled in R.

### How it works

1. **Initial BLASTp**: query protein vs. subject proteome.
2. **Filtering**: identity ≥ 30%, coverage ≥ 21%, E-value < 1e-10.
3. **Reiterative search**: each filtered hit becomes a new query.
4. **Deduplication**: tracks used queries to avoid redundant searches.
5. **Excel output**: one sheet per query protein + a final deduplicated sheet.

### Requirements
```r
install.packages(c("tidyverse", "openxlsx"))
BiocManager::install(c("Biostrings", "ape"))
```
NCBI BLAST+ must be installed and `blastp` accessible.

### Usage
Edit the configuration block at the top of the script:
```r
query_fasta_file   <- "your_query.fasta"
subject_fasta_file <- "your_subject.faa"
blastp_executable  <- "blastp"
output_excel       <- "Blast_prot_local.xlsx"
```
Then run:
```r
source("reiterative_BLASTp_local.R")
```

### Output
An Excel file with one sheet per query protein searched, plus a `Unique_results` sheet with the final deduplicated hits.

---

## Step 2 — `extract_unique_blast_sequences.R`

Collects `.xlsx` BLAST results from a multi-species results directory, deduplicates hits per species, extracts the matching sequences from each proteome database, and merges everything into a single annotated FASTA.

### How it works

1. **Scans** a directory where each subfolder contains `.xlsx` BLAST results for one species.
2. **Deduplicates** by `Subject` column, saving `datos_unicos.xlsx` per folder.
3. **Extracts** matching sequences from each species' protein FASTA.
4. **Annotates** each sequence header with the species name.
5. **Merges** all sequences into a single output FASTA.

### Requirements
```r
install.packages(c("readxl", "writexl", "dplyr"))
BiocManager::install("Biostrings")
```

### Expected directory structure
```
results/reiterative_blastp_results/
├── Species_A/
│   ├── blast_round1.xlsx
│   └── blast_round2.xlsx
├── Species_B/
│   └── blast_round1.xlsx
```

### Usage
Edit the configuration block at the top of the script:
```r
results_dir <- "results/reiterative_blastp_results"

species_db <- list(
  "Species_A" = "data/raw/species_a_protein.faa",
  "Species_B" = "data/raw/species_b_protein.faa"
)

query_fasta  <- "data/raw/my_query_proteins.fasta"  # optional
output_fasta <- "todas_las_secuencias.faa"
```
Then run:
```r
source("extract_unique_blast_sequences.R")
```

### Output
| File | Description |
|------|-------------|
| `datos_unicos.xlsx` | Deduplicated hits per species folder |
| `todas_las_secuencias.faa` | Final merged FASTA, species-annotated |

---

## 🔬 Biological context

Developed to identify and consolidate **protein family members** across multiple species proteomes. The reiterative approach captures distant homologs missed by a single BLAST run, while the extraction step prepares sequences for **multiple sequence alignment and phylogenetic tree construction**.

---

## 👤 Author

**Sergio Civera Arroyo**
Biologist & Bioinformatician | MSc Bioinformatics
civeraroysergio@gmail.com
