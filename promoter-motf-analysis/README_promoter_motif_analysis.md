# Promoter Cis-Regulatory Motif Analysis

A statistical pipeline in R for detecting and characterizing **cis-regulatory motifs** in gene promoter regions. The script identifies whether specific sequence motifs appear more frequently than expected by chance, using **Monte Carlo simulation** to build a null distribution.

---

## What does this script do?

1. **Reads input data**: Parses a FASTA file with promoter sequences and a CSV/TSV file defining the motifs to search.
2. **Motif search**: Scans each promoter for all motif variants, with optional IUPAC degenerate base support and reverse-complement search.
3. **Monte Carlo simulation**: Generates a null distribution (1000+ iterations) using one of three null models: mononucleotide frequencies, dinucleotide frequencies, or sequence permutation.
4. **Statistical enrichment**: Calculates enrichment scores and empirical p-values for each motif, corrected for multiple testing (Benjamini-Hochberg or Bonferroni).
5. **Positional analysis**: Computes absolute and relative positional distributions of motif occurrences across promoters.
6. **Outputs**: Exports summary tables (CSV) and publication-ready plots (PNG).

---

## Requirements

### R packages
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "stringr",
                   "readr", "purrr", "scales", "ggridges"))
```

> No Bioconductor packages required — all dependencies are on CRAN.

---

## Input files

### 1. Promoter FASTA file
Standard FASTA format (`.fa`, `.fasta`, `.txt`). One sequence per promoter.

```
>GeneID_1
ATCGATCGATCGATCG...
>GeneID_2
GCTAGCTAGCTAGCTA...
```

- Headers must contain a unique identifier after `>`
- Sequences can span multiple lines
- IUPAC ambiguous bases (N, R, Y...) are accepted and documented

### 2. Motif definition file
CSV or TSV file (auto-detected) with at least two columns:

| Column | Description |
|--------|-------------|
| `motif_id` | Motif name or family (e.g. `ABRE`, `DRE`, `GBOX`) |
| `sequence` | Exact motif sequence (e.g. `ACGTG`, `CACGTG`) |

Optional columns: `group`, `function`, `source`.

```
motif_id,sequence
ABRE,ACGTG
ABRE,CACGTG
DRE,GCCGAC
GBOX,CACGTG
```

> Multiple rows with the same `motif_id` represent equivalent variants of the same regulatory element.

---

## Usage

1. Edit the configuration block at the top of the script:

```r
# ── USER CONFIGURATION ──────────────────────────────────────────────────────
ruta_fasta   <- "promotores.fa"       # FASTA file with promoter sequences
ruta_motivos <- "motivos.tsv"         # CSV/TSV file with motif definitions
ruta_salida  <- "resultados_motivos"  # Output folder (created automatically)

n_sim         <- 2500     # Monte Carlo iterations (minimum 1000 recommended)
semilla       <- 42       # Random seed for reproducibility

usar_iupac      <- TRUE   # Expand IUPAC codes in motifs
buscar_rev_comp <- TRUE   # Search reverse-complement strand too

modelo_nulo   <- "dinucleotidica"  # Null model: "mononucleotidica", "dinucleotidica", "permutacion"
metodo_ajuste <- "BH"              # P-value correction: "BH" (FDR) or "bonferroni"
umbral_padj   <- 0.05              # Significance threshold
# ────────────────────────────────────────────────────────────────────────────
```

2. Set your working directory to the folder containing the input files.

3. Run the script:
```r
source("analisis_motivos_promotores.R")
```

---

## Output

All results are saved in the output folder (`ruta_salida`):

| File | Description |
|------|-------------|
| `tabla_resumen_motivos.csv` | Main results table — one row per motif with enrichment and p-values |
| `tabla_ocurrencias_reales.csv` | One row per motif occurrence found |
| `tabla_posiciones_relativas.csv` | Relative positions (0–1) of each occurrence |
| `tabla_bins_posicionales_relativos.csv` | Binned relative positional distribution |
| `tabla_bins_posicionales_absolutos.csv` | Binned absolute positional distribution (bp) |
| `resumen_analisis.txt` | Human-readable summary of parameters and significant results |
| `*.png` | Enrichment plots and positional distribution charts |

The `resultados` object returned by the pipeline can also be explored directly in R:
```r
resultados$resumen_final     # Main statistics per motif
resultados$ocurrencias       # Full occurrences table
resultados$posiciones        # Positional data
resultados$simulaciones      # Null distributions per motif
```

---

## Statistical approach

- **Enrichment score**: Ratio of observed vs. mean expected occurrences under the null model.
- **Empirical p-value**: Fraction of simulations where expected ≥ observed.
- **Two metrics tested**: number of promoters containing the motif, and total number of occurrences.
- **Multiple testing correction**: Applied independently to both metrics.

---

## Biological context

Developed to characterize **stress-responsive cis-regulatory elements** in the promoters of genes of interest in plant species. Useful for any organism where promoter sequences and candidate regulatory motifs are available.

---

## Author

**Sergio Civera Arroyo**  
Biologist & Bioinformatician

