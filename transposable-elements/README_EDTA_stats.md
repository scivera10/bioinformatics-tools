# EDTA TE Annotation Report Generator

An R script that parses a **GFF3 file produced by EDTA/panEDTA** and generates a structured report of transposable element (TE) annotations, replicating the style of the original RepeatMasker summary report.

\---

## What does this script do?

Given an EDTA-annotated genome GFF3 file and the total genome size, the script produces three tables:

1. **Repeat Classes** — TE counts and masked base pairs grouped by class (LTR, TIR, LINE, SINE, etc.) and subclass (Gypsy, Copia, Mutator, etc.), including % of genome masked.
2. **Repeat Stats** — Counts and masked bp broken down by TE family name.
3. **By Sequence** — Per-chromosome/scaffold summary of TE count and masked bp.

All three tables are printed to the console and exported as CSV files.

\---

## Requirements

### R packages

```r
install.packages(c("data.table", "dplyr", "tidyr", "stringr", "beepr"))
```

### Input

* A **GFF3 file** produced by [EDTA](https://github.com/oushujun/EDTA) or panEDTA
* The **total genome size in base pairs**

\---

## Usage

### From the command line (recommended)

```bash
Rscript EDTA\_stats.R input.gff3 715869006 output\_prefix
```

|Argument|Description|
|-|-|
|`input.gff3`|Path to the EDTA GFF3 annotation file|
|`715869006`|Total genome size in bp|
|`output\_prefix`|Prefix for output CSV files (default: `TE\_report`)|

### From RStudio

Edit the fallback values at the top of the configuration block:

```r
gff\_file    <- "input.gff3"    # your GFF3 file
genome\_size <- 1000000000      # your genome size in bp
out\_prefix  <- "TE\_report"
```

Then run the script.

\---

## Output

|File|Description|
|-|-|
|`<prefix>\_repeat\_classes.csv`|TE counts by class and superclass|
|`<prefix>\_repeat\_stats.csv`|TE counts by family name|
|`<prefix>\_by\_sequence.csv`|TE counts per chromosome/scaffold|

\---

## Methodological note

> The original EDTA report is generated from the RepeatMasker `.out` file, which resolves overlapping annotations (highest score wins). The GFF3 retains \*\*all\*\* annotations including overlapping ones. As a result, counts from this script may be slightly higher (\~5%) than the original report. Base pair counts for homology-based elements are exact; structural LTR retrotransposons use the element body length (excluding TSDs).

\---

## Supported TE classes

|GFF3 type|Report class|
|-|-|
|`Gypsy\_LTR\_retrotransposon`|LTR/Gypsy|
|`Copia\_LTR\_retrotransposon`|LTR/Copia|
|`LTR\_retrotransposon`|LTR/unknown|
|`L1\_LINE\_retrotransposon`|LINE/L1|
|`RTE\_LINE\_retrotransposon`|LINE/RTE|
|`tRNA\_SINE\_retrotransposon`|SINE/tRNA|
|`CACTA\_TIR\_transposon`|TIR/CACTA|
|`Mutator\_TIR\_transposon`|TIR/Mutator|
|`PIF\_Harbinger\_TIR\_transposon`|TIR/PIF\_Harbinger|
|`Tc1\_Mariner\_TIR\_transposon`|TIR/Tc1\_Mariner|
|`hAT\_TIR\_transposon`|TIR/hAT|
|`helitron`|nonTIR/helitron|
|`DIRS\_YR\_retrotransposon`|nonLTR/DIRS\_YR|
|`repeat\_fragment`|repeat\_fragment|
|`low\_complexity`|low\_complexity|

\---

## Biological context

Developed for the analysis of **transposable element landscapes** in plant genomes annotated with EDTA. Useful for any eukaryotic genome where EDTA has been used for TE annotation.

\---

## Author

**Sergio Civera Arroyo**  
Biologist \& Bioinformatician


