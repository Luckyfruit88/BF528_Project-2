# BF528 Project 2: RNA-seq Pipeline and Differential Expression Analysis

## 📌 Project Overview
This repository contains an automated RNA-seq pipeline and downstream statistical analysis for BU BF528 (Genomic Data Analysis). The workflow processes raw paired-end reads, aligns them to the human reference genome, quantifies gene expression, and supports differential expression analysis.

The analysis goal is to process data from GSE190725 and reproduce key findings from the Nature Communications study: *"The type 1 diabetes gene TYK2 regulates β-cell development and its responses to interferon-α"*.

## 🧬 Biological Context
- **Condition:** Human induced pluripotent stem cells (hiPSCs), Stage 5 endocrine progenitors.
- **Experimental Group:** TYK2 Knockout (KO)
- **Control Group:** Wild Type (WT)
- **Objective:** Identify transcriptomic changes caused by TYK2 ablation.

---

## 🛠️ Pipeline Architecture
Upstream processing is implemented in **Nextflow DSL2** and containerized.

### 1. Data Preprocessing & Alignment (Nextflow)
- **Quality Control:** `FastQC` (container: `ghcr.io/bf528/fastqc:latest`)
- **Genome Indexing & Alignment:** `STAR` (container: `ghcr.io/bf528/star:latest`) against GRCh38 primary assembly
- **Quantification:** `VERSE` (container: `ghcr.io/bf528/verse:latest`) with `gencode.v45.primary_assembly.annotation.gtf`
- **Matrix Aggregation:** Python + `pandas` merge of per-sample `.exon.txt` outputs
- **QC Aggregation:** `MultiQC` (container: `ghcr.io/bf528/multiqc:latest`)

### 2. Downstream Statistical Analysis (R / Rmarkdown)
- **Filtering:** Remove low-count genes
- **Normalization & QC:** `rlog` / `vst`, PCA, sample distance heatmaps
- **Differential Expression:** `DESeq2` (KO vs WT)
- **Pathway Enrichment:** `fgsea` with MSigDB C2

---

## ⚙️ Runtime Parameters
Values are defined in `nextflow.config`:
- `params.reads`: `/projectnb/bf528/materials/project-2-rnaseq/full_files/*_R{1,2}.fastq.gz`
- `params.genome`: `/projectnb/bf528/materials/project-2-rnaseq/refs/GRCh38.primary_assembly.genome.fa`
- `params.gtf`: `/projectnb/bf528/materials/project-2-rnaseq/refs/gencode.v45.primary_assembly.annotation.gtf`
- `params.outdir`: `$projectDir/results/`

## 📂 Repository Structure
```text
├── nextflow.config
├── main.nf
├── project-2-report.Rmd
├── project-2-report.html
├── gitignore.txt
└── README.md
```

## ❗ Workflow Dependency Note
`main.nf` currently expects helper scripts at:
- `bin/gtf_to_symbol.py`
- `bin/merge_counts.py`

Create/populate these files before running the Nextflow workflow.