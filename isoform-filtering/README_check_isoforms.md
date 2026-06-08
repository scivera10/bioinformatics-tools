# Isoform Filtering — One Representative Protein per Locus

An R function that selects a single **representative protein isoform per gene locus** from a GFF3 annotation and its associated protein FASTA. Designed for multi-species comparative genomics workflows where redundant isoforms need to be removed before downstream analyses (phylogenetics, BLAST, etc.).

\---

## What does this script do?

Given a GFF3 annotation and a protein FASTA file, for each gene locus it selects the single best isoform using a three-tier priority criterion:

|Priority|Criterion|Rationale|
|-|-|-|
|1st|**Longest protein (aa)**|Direct proxy of CDS completeness|
|2nd|**Longest total CDS (bp)**|Tiebreaker summing all CDS features|
|3rd|**Lowest isoform number**|Reproducibility when all else is equal|

The script is compatible with NCBI GFF3 formats from different species, automatically detecting whether to use `locus\_tag` or `gene` as the locus identifier.

\---

## Requirements

### R packages

```r
install.packages(c("dplyr", "cli", "beepr"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))
```

\---

## Input files

|File|Description|
|-|-|
|`species.gff`|GFF3 annotation file (NCBI format recommended)|
|`species\_protein.faa`|Protein FASTA file matching the GFF3 annotation|

The script automatically handles GFF3 files with or without `locus\_tag` column (e.g. chromosome-level vs scaffold-level NCBI assemblies).

\---

## Usage

1. Source the script to load the `check\_isoforms()` function:

```r
source("check\_isoforms\_v4.R")
```

2. Call the function with your input files:

```r
result <- check\_isoforms(
  gff\_file   = "species.gff",
  faa\_file   = "species\_protein.faa",
  output\_dir = "output/species"
)
```

3. Or edit the execution blocks at the bottom of the script and toggle `FALSE -> TRUE`.

\---

## Output

|File|Description|
|-|-|
|`<input>\_cleaned\_isoforms.faa`|Filtered FASTA with one protein per locus|

The function returns an R list:

```r
result$table       # Selection table with chosen protein\_id per locus
result$fasta       # Filtered AAStringSet object
result$output\_faa  # Path to the output FASTA
result$stats       # Summary statistics tibble
```

### Summary statistics

|Metric|Description|
|-|-|
|`total\_mrna\_features`|Total mRNA features in the GFF3|
|`total\_loci`|Unique gene loci detected|
|`loci\_with\_protein`|Loci with a valid protein selected|
|`loci\_lost`|Loci with no valid protein in the FASTA|
|`fasta\_input\_seqs`|Sequences in input FASTA|
|`fasta\_output\_seqs`|Sequences in output FASTA|
|`isoforms\_removed`|Isoforms discarded|

\---

## Biological context

Developed for **comparative genomics and phylogenetic analyses** across multiple plant species. Redundant isoforms inflate BLAST hit counts and distort phylogenetic trees — this script ensures each gene is represented exactly once before running inter-species analyses.

Compatible with NCBI RefSeq protein datasets for any eukaryotic species.

\---

## Author

**Sergio Civera Arroyo**
Biologist \& Bioinformatician

