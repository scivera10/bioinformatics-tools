# Sanger Validation Pipeline

A two-step Bash pipeline for validating NGS variants by extracting reads from a genomic region of interest and realigning them against a Sanger consensus sequence. Useful for confirming variants identified in whole-genome or targeted sequencing experiments.

\---

## Pipeline overview

```
BAMs from NGS mapping
        │
        ▼
┌──────────────────────────────────┐
│  extraccion\_region\_interes.sh    │  ── Step 1: extract reads from target region
└──────────────────────────────────┘
        │
        ▼
  Per-sample FASTQs (R1 + R2)
        │
        ▼
┌──────────────────────────────────┐
│  alineamiento\_sanger.sh          │  ── Step 2: align against Sanger reference
└──────────────────────────────────┘
        │
        ▼
  Per-sample sorted + indexed BAMs (ready for IGV or variant calling)
```

\---

## Requirements

### External tools

* [samtools](http://www.htslib.org/) ≥ 1.13
* [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) ≥ 2.4

Both must be installed and accessible from the `PATH`.

### Bowtie2 index

The Sanger reference must be indexed before running Step 2:

```bash
bowtie2-build sanger\_reference.fasta sanger\_referencia/sanger\_index
```

\---

## Expected directory structure

```
project\_root/
├── Sample1/
│   └── resecuenciacion/
│       └── 2\_bowtie2\_mapping/
│           └── Sample1.sorted.bam       ← input BAM (indexed)
├── Sample2/
│   └── ...
├── sanger\_referencia/
│   └── sanger\_index.\*                   ← Bowtie2 index files
└── NGS\_SANGER/                          ← created automatically by Step 1
    ├── Sample1\_SANGER/
    │   ├── Sample1\_R1.fq.gz
    │   ├── Sample1\_R2.fq.gz
    │   └── ...
    └── Sample2\_SANGER/
        └── ...
```

\---

## Step 1 — `extraccion\_region\_interes.sh`

Extracts paired-end reads from a specific genomic region from each sample's BAM file and converts them to FASTQ.

### Configuration

Edit the configuration block at the top of the script:

```bash
REGION="chrX:1000000-1005000"   # Target region in samtools format (chr:start-end)

MUESTRAS=(
    Sample1
    Sample2
    Sample3
)
```

### Usage

```bash
cd /path/to/project\_root
bash extraccion\_region\_interes.sh
```

### Output per sample (in `NGS\_SANGER/<SAMPLE>\_SANGER/`)

|File|Description|
|-|-|
|`<SAMPLE>\_R1.fq.gz`|Forward reads from the target region|
|`<SAMPLE>\_R2.fq.gz`|Reverse reads from the target region|
|`singletons.fastq.gz`|Unpaired reads|
|`flagstat.txt`|Alignment statistics for the extracted region|

\---

## Step 2 — `alineamiento\_sanger.sh`

Aligns the extracted FASTQ reads against the Sanger consensus reference using Bowtie2. Produces a sorted, indexed BAM and an alignment log per sample.

### Configuration

Edit the configuration block at the top of the script:

```bash
DIRECTORIO\_BASE="."                        # Project root directory
DIR\_REFERENCIA="${DIRECTORIO\_BASE}/sanger\_referencia"
INDICE\_SANGER="${DIR\_REFERENCIA}/sanger\_index"
THREADS=4

MUESTRAS=(
    Sample1
    Sample2
    Sample3
)
```

### Usage

```bash
cd /path/to/project\_root
bash alineamiento\_sanger.sh
```

### Output per sample (in `NGS\_SANGER/<SAMPLE>\_SANGER/`)

|File|Description|
|-|-|
|`<SAMPLE>.sanger.sorted.bam`|Sorted BAM aligned to Sanger reference|
|`<SAMPLE>.sanger.sorted.bam.bai`|BAM index (for IGV or variant calling)|
|`<SAMPLE>.sanger.bowtie2.log`|Bowtie2 alignment statistics|

\---

## Biological context

Developed for validating variants identified by NGS in targeted genomic regions. After alignment against a Sanger-validated consensus, the resulting BAMs can be inspected in **IGV** or used for local variant calling to confirm or reject candidate variants.

\---

## Author

**Sergio Civera Arroyo**
Biologist \& Bioinformatician

