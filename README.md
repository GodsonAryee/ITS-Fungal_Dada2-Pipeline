````
# ITS Fungal DADA2 Pipeline (R)

A reproducible **ITS (fungal) amplicon sequencing pipeline** using **DADA2**, **cutadapt**, and **phyloseq** for processing paired-end Illumina data.  
This workflow performs **quality filtering, primer removal, denoising, chimera removal, taxonomic assignment (UNITE), and downstream ecological analysis**, producing publication-ready outputs.

Repository script: `8.7.2025_ITS.R`

---

## Overview

This pipeline is designed for **paired-end ITS1 fungal sequencing data** and includes:

1. Pre-filtering reads to remove ambiguous bases  
2. Primer detection and removal using **cutadapt**  
3. DADA2 denoising (ASV inference)  
4. Chimera removal  
5. Taxonomic assignment using the **UNITE fungal reference database**  
6. Construction of **phyloseq objects**  
7. Read tracking and quality summaries  
8. Genus-level relative abundance visualization  

---

## Requirements

### Software
- **R ≥ 4.2**
- **cutadapt ≥ 4.0**
- UNIX/Linux or macOS environment recommended

### R packages
- `dada2`
- `ShortRead`
- `Biostrings`
- `phyloseq`
- `ggplot2`
- `dplyr`
- `data.table`

Install core dependencies:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("dada2", "ShortRead", "Biostrings", "phyloseq"))
install.packages(c("ggplot2", "dplyr", "data.table"))
````

---

## Input Data Requirements

### FASTQ files

Paired-end reads must follow this naming convention:

```
SampleID (paired)_ITS1_R1.fastq.gz
SampleID (paired)_ITS1_R2.fastq.gz
```

The pipeline automatically:

* Detects forward (`_R1`) and reverse (`_R2`) reads
* Extracts sample names
* Ensures correct pairing

### Metadata

A tab-delimited metadata file is required:

```
METADATA.txt
```

* Row names must match sample IDs
* Used to construct the phyloseq object

---

## Primer Information

The pipeline removes ITS primers using **cutadapt**:

* **Forward (ITS1F):** `CTTGGTCATTTAGAGGAAGTAA`
* **Reverse (ITS2):** `GCTGCGTTCTTCATCGATGC`

All primer orientations are checked prior to trimming to confirm successful removal.

---

## How to Run

1. Edit paths in the script:

```r
path <- "/path/to/ITS_fastq_files/"
cutadapt <- "/path/to/cutadapt"
unite.ref <- "/path/to/UNITE_reference.fasta"
```

2. Run the script:

```bash
Rscript 8.7.2025_ITS.R
```

---

## Pipeline Steps (Detailed)

### 1. Pre-filtering

* Removes reads with Ns
* Filters low-quality reads
* Minimum length: **50 bp**

### 2. Primer Removal

* Uses **cutadapt** with paired-end synchronization
* Verifies primer removal before and after trimming

### 3. DADA2 Processing

* Error rate learning
* Dereplication
* ASV inference
* Paired-read merging
* Chimera removal (consensus method)

### 4. Read Tracking

Tracks reads through each stage:

* Input
* Filtered
* Denoised (forward/reverse)
* Merged
* Non-chimeric

Saved as:

```
read_tracking.csv
```

### 5. Taxonomic Assignment

* UNITE fungal reference database
* Reverse-complement matching enabled
* Outputs taxonomy table

### 6. Phyloseq Object Construction

* ASVs relabeled as `ASV1`, `ASV2`, ...
* Non-target taxa removed:

  * Eukaryota
  * Chloroplast
  * Mitochondria
* Saves cleaned phyloseq object

---

## Outputs

| File                         | Description                      |
| ---------------------------- | -------------------------------- |
| `seqtab_nochim.rds`          | ASV table (non-chimeric)         |
| `taxa_assignments.rds`       | Taxonomy assignments             |
| `read_tracking.csv`          | Read counts through pipeline     |
| `ps_water.rds`               | Clean phyloseq object            |
| `Reads_per_sample_Water.csv` | Read depth per sample            |
| Genus-level bar plots        | Relative abundance visualization |

---

## Visualization

The pipeline generates **genus-level relative abundance bar plots** using `ggplot2`, with:

* Relative abundance normalization
* Sample-wise stacking
* Rotated x-axis labels for clarity

---

## Notes & Best Practices

* ITS regions are highly variable in length; avoid aggressive truncation
* Always inspect quality profiles before adjusting filters
* UNITE reference version should match your study year
* For large datasets, HPC execution is recommended

---

## Citation

If you use this pipeline, please cite:

Callahan BJ et al. (2016).
**DADA2: High-resolution sample inference from Illumina amplicon data.**
*Nature Methods*, 13:581–583.
---

## Author

Godson Aryee

* GitHub: `https://github.com/GodsonAryee`
* Email: `godson.aryee@ndsu.edu`
