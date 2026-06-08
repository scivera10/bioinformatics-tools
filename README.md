# Bioinformatics Tools — Sergio Civera Arroyo

A collection of R scripts for genomic and bioinformatic analyses, developed during research projects in plant genomics and comparative genomics. Each folder contains a self-contained pipeline or tool with its own documentation.

\---

## Repository structure

```
bioinformatics-tools/
├── blast/                        # Reiterative BLASTp pipeline
├── isoform-filtering/            # Representative isoform selection
├── promoter-extraction/          # Promoter sequence extraction
├── promoter-motif-analysis/      # Cis-regulatory motif enrichment
├── transposable-elements/        # EDTA TE annotation report
└── utils/                        # General-purpose utilities
```

\---

## Tools

### [`blast/`](./blast/)

**Reiterative BLASTp Pipeline** — Two-step pipeline for identifying protein homologs across multiple species proteomes.

* `reiterative\_BLASTp\_local.R` — Runs iterative local BLASTp searches using each hit as a new query, expanding the search until no new homologs are found. Results exported to Excel.
* `extract\_unique\_blast\_sequences.R` — Collects and deduplicates hits from a multi-species BLAST results directory, extracts matching sequences from each proteome, and merges them into a single species-annotated FASTA.

\---

### [`isoform-filtering/`](./isoform-filtering/)

**Isoform Filtering** — Selects one representative protein per gene locus from a GFF3 annotation and its associated protein FASTA, using a three-tier criterion: longest protein → longest CDS → lowest isoform number. Compatible with NCBI GFF3 formats across species.

\---

### [`promoter-extraction/`](./promoter-extraction/)

**Promoter Sequence Extraction** — Extracts promoter regions from a reference genome including the 5'UTR but trimming any overlap with the CDS, on both strands. Also exports GFF3 and BED files for visualization in IGV. Output feeds directly into `promoter-motif-analysis/`.

\---

### [`promoter-motif-analysis/`](./promoter-motif-analysis/)

**Cis-Regulatory Motif Enrichment** — Statistical pipeline for detecting enriched sequence motifs in promoter regions. Uses Monte Carlo simulation (1000+ iterations) to build a null distribution and computes empirical p-values corrected for multiple testing. Supports IUPAC degenerate bases, reverse-complement search, and three null models (mononucleotide, dinucleotide, permutation).

\---

### [`transposable-elements/`](./transposable-elements/)

**EDTA TE Annotation Report** — Parses a GFF3 file produced by EDTA/panEDTA and generates a RepeatMasker-style summary report with TE counts, masked base pairs, and genome coverage — broken down by class, family, and chromosome.

\---

### [`utils/`](./utils/)

General-purpose utilities:

* `separar\_seq\_fasta.R` — Splits a multi-sequence FASTA into individual files, one per sequence, named after the sequence header.
* `extract\_selected\_proteins.R` — Extracts a subset of protein sequences from a FASTA database based on accession IDs provided in an Excel file.

\---

## General requirements

All scripts are written in **R**. Package requirements are listed in each tool's README. Most tools use a combination of:

* [Bioconductor](https://bioconductor.org/) packages: `Biostrings`, `GenomicRanges`, `rtracklayer`, `Rsamtools`
* CRAN packages: `tidyverse`, `dplyr`, `ggplot2`, `openxlsx`, `readxl`
* External tools: [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (for `blast/`)

\---

## Author

**Sergio Civera Arroyo**  
Biologist \& Bioinformatician

