# Promoter Sequence Extraction (5'UTR included, CDS excluded)

An R script for extracting **promoter sequences** from a reference genome, including the 5'UTR region but precisely trimming any overlap with the CDS. Also generates GFF3 and BED files ready for visualization in **IGV**.

\---

## What does this script do?

1. **Reads a gene list** from an Excel file containing Entrez gene IDs.
2. **Imports a genome annotation** (GFF/GFF3) and filters the genes of interest.
3. **Defines promoter regions** as N bp upstream of the TSS (default: 1500 bp).
4. **Trims CDS overlap**: adjusts promoter coordinates so they include the 5'UTR but stop exactly at the start of the coding sequence — on both + and − strands.
5. **Validates** the final regions (removes invalid ranges, checks for residual CDS overlap).
6. **Extracts sequences** from the reference genome FASTA.
7. **Exports IGV-ready files**: GFF3 and BED files for genes, promoters, and CDS.

\---

## Requirements

### R packages

```r
install.packages(c("tidyr", "dplyr", "readxl"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer", "Rsamtools", "Biostrings"))
```

### External requirements

* Reference genome FASTA **indexed** with `samtools faidx`:

```bash
samtools faidx genome.fna
```

\---

## Input files

|File|Description|
|-|-|
|`gene\_table.xlsx`|Excel file with a column named `Entrez ID` containing gene identifiers|
|`annotation.gff`|GFF/GFF3 genome annotation file|
|`gr\_gff\_filtrado\_CG8.gff3`|Manually curated GFF3 (optional — for adding genes missing from the annotation)|
|`genome.fna`|Reference genome FASTA (must be indexed with `samtools faidx`)|

\---

## Usage

1. Edit the configuration block at the top of the script:

```r
# ── USER CONFIGURATION ──────────────────────────────────────────────────────
gene\_table\_xlsx   <- "gene\_table.xlsx"
gff\_file          <- "annotation.gff"
gff\_filtered\_CG8  <- "gr\_gff\_filtrado\_CG8.gff3"
genome\_fasta      <- "genome.fna"
upstream\_bp       <- 1500    # bp upstream of TSS
# ────────────────────────────────────────────────────────────────────────────
```

2. Set your working directory to the folder containing all input files.
3. Run the script:

```r
source("Promotores\_sin\_CDS.R")
```

> \*\*Note\*\*: After step 3, the script exports a `gr\_gff\_filtrado.gff3` file for manual review. If you need to add genes not found in the annotation, add them manually, save as `gr\_gff\_filtrado\_CG8.gff3`, and continue running the rest of the script.

\---

## Output

|File|Description|
|-|-|
|`secuencias\_promotores\_sin\_CDS\_con\_UTR.fa`|FASTA file with final promoter sequences|
|`gr\_gff\_filtrado.gff3`|Intermediate filtered GFF3 for manual review|
|`IGV\_genes\_oleosinas.gff3`|GFF3 of target genes for IGV|
|`IGV\_promotores\_sin\_CDS\_con\_5UTR.gff3`|GFF3 of final promoter regions for IGV|
|`IGV\_CDS\_oleosinas.gff3`|GFF3 of CDS regions for IGV|
|`IGV\_genes\_oleosinas.bed`|BED of target genes for IGV|
|`IGV\_promotores\_sin\_CDS\_con\_5UTR.bed`|BED of final promoter regions for IGV|
|`IGV\_CDS\_oleosinas.bed`|BED of CDS regions for IGV|

\---

## Biological context

Developed to extract promoter sequences for downstream **cis-regulatory motif analysis** (see `promoter-motif-analysis/`). The precise exclusion of the CDS while retaining the 5'UTR ensures that regulatory elements near the transcription start site are captured without including coding sequence.

Compatible with any annotated genome where gene models include CDS features in the GFF.

\---

## Author

**Sergio Civera Arroyo**  
Biologist \& Bioinformatician

