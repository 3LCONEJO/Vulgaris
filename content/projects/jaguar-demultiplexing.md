---
title: "JAGUAR-DEMULTIPLEXING"
date: 2026-06-26
draft: false
status: "complete"
summary: "Given the information in my garden, an algorithm helps me infer which berry belongs to which bush. In our data, it tells us which cell belongs to which individual."
tags: ["scATAC-seq", "JAGUAR", "Demultiplexing", "Bash", "LSF"]
github: "https://github.com/JAGUAR-LATAM/Jaguar-Fastdemux/"
stack: ["popscle", "bcftools", "samtools", "bash"]
---

> DIEGO RAMIREZ HELP ME WITH THIS SCRIPT

> THE GITHUB SENDO YOU TO THE FASTDEMUX VERSION OF THE PIPELINE BUT ESSENTIALY IS THE SAME, THE ONLY THING THAT CHANGE IS THE FASTEDMUX PART AND REPLACE IT WITH DEMUXLET

## The Problem

Think of our dataset like a massive community garden. If we want to understand the differences in immunity across populations, we need to know exactly where every single piece of produce came from. 

Since Diego and I are handling the ATAC data, we must precisely label each cell with its individual sample identity. Without this fundamental step, we cannot perform accurate differential accessibility analyses downstream.

---

## Approach: Demuxlet

After several trials evaluating different architectures, we concluded that `popscle demuxlet` is the most robust tool for this specific data. It leverages genotype likelihoods to map individual cells back to their donors of origin, identifying both singlets and doublets.

To ensure the script runs efficiently and doesn't crash on the cluster, the pipeline heavily filters and reorganizes the genetic data *before* feeding it to the demultiplexer.

### Pipeline Logic

We use a strict step-by-step logic to reduce the computational load. By restricting our search space only to ATAC peaks and strictly biallelic SNPs, we save massive amounts of memory.

`[Raw BAM & Joint VCF]` ->

`[Subset to Donors]` ->

`[Filter Biallelic SNPs]` ->

`[Restrict to ATAC Peaks]` ->

`[Reorder Contig Headers]` ->

`[Demuxlet Assignment]`

---

## Implementation Guide

Below is the complete, idempotent bash script designed for an LSF cluster environment. Because every step checks if its output already exists, you can safely restart the job if a node fails without losing your progress.

### 1. Cluster Setup and Configuration
We start by defining the LSF job parameters, loading the necessary modules (`bcftools`, `samtools`, and the Singularity-based `demuxafy`), and setting up our global paths.

```bash
#!/bin/bash
#BSUB -G jaguar_analysis
#BSUB -J demuxlet_scATAC
#BSUB -q normal
#BSUB -n 16
#BSUB -R "select[mem>64000] rusage[mem=64000] span[hosts=1]"
#BSUB -M 64000
#BSUB -q long
#BSUB -o logs/demuxlet_%J.out
#BSUB -e logs/demuxlet_%J.err
#BSUB -W 48:00

set -euo pipefail

echo "--- Loading software modules ---"
module load cellgen/bcftools/1.23
module load HGI/common/demuxafy/2.1.0 

THREADS="${LSB_DJOB_NUMPROC:-16}"
CONFIG_TSV="config_pipeline.tsv"
GLOBAL_ANALYSIS_DIR="/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jv8/Demuxlet-scATAC-bge09"
JOINT_DONOR_VCF="/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/ks37/gt/qc/BGE09/filtered/output/bge09_exome_biallelic_pass_only.snp_indel_le50.setid.autosomes_only.JAGID.dedup.vcf.gz"
BED_FOLDER="/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jv8/beds"

mkdir -p logs
log() { echo "[$(date '+\%Y-\%m-\%d \%H:\%M:\%S')]$*"; }
```

### 2. Initialization and BAM Verification

The script loops through the configuration file. For each sample, it verifies that the BAM file exists and generates an index if one is missing. It then extracts the chromosome order directly from the BAM header, which is strictly required by popscle later on.

```R
tail -n +2 "$CONFIG_TSV" \vert{} while IFS=$'\t' read -r SAMPLE_ID BAM_FILE BARCODE_FILE DONOR_LIST_FILE NUM_DONORS
do
  BASE_DIR="${GLOBAL_ANALYSIS_DIR}/${SAMPLE_ID}"
  DEMUXLET_DIR="${BASE_DIR}/demuxlet_out"
  DONOR_VCF_DIR="${BASE_DIR}/donor_vcf_prep"
  mkdir -p "$DEMUXLET_DIR" "$DONOR_VCF_DIR"

  # STEP -1: Verify BAM exists and is indexed
  if [ ! -f "${BAM_FILE}.bai" ]; then
    log "BAM index missing — creating it..."
    samtools index -@ "$THREADS" "$BAM_FILE"
  fi

  # STEP 0: Extract chromosome order from BAM header (@SQ lines)
  BAM_ORDER_FILE="${DONOR_VCF_DIR}/bam.chrom.order.txt"
  if [ ! -f "$BAM_ORDER_FILE" ]; then
    log "Extracting chromosome order from BAM header..."
    samtools view -H "$BAM_FILE" \vert{} grep '^@SQ' \vert{} cut -f2 \vert{} sed 's/SN://' > "$BAM_ORDER_FILE"
  fi
```

### 3. VCF Filtering (Donors & Biallelic SNPs)

Running pileups on whole-genome VCFs is highly inefficient. Here, we subset the massive joint VCF specifically to the donors in our pool, and then aggressively filter out indels and multi-allelic sites.

```R
# STEP 1: Subset joint germline VCF to the donors of interest
  SUBSET_VCF="${DONOR_VCF_DIR}/donors.subset.vcf.gz"
  if [ ! -f "$SUBSET_VCF" ]; then
    log "Subsetting VCF to donors in $DONOR_LIST_FILE ..."
    bcftools view --threads "$THREADS" -S "$DONOR_LIST_FILE" "$JOINT_DONOR_VCF" -Oz -o "${DONOR_VCF_DIR}/tmp.subset.vcf.gz"
    bcftools sort "${DONOR_VCF_DIR}/tmp.subset.vcf.gz" -Oz -o "$SUBSET_VCF"
    tabix -p vcf "$SUBSET_VCF"
  fi

  # STEP 1.5: Filter to strictly biallelic SNPs
  BIALLELIC_SNPS_VCF="${DONOR_VCF_DIR}/donors.subset.biallelic_snps.vcf.gz"
  if [ ! -f "$BIALLELIC_SNPS_VCF" ]; then
    log "Filtering to biallelic SNPs..."
    bcftools view --threads "$THREADS" -v snps -m2 -M2 "$SUBSET_VCF" -Oz -o "${DONOR_VCF_DIR}/tmp.biallelic_snps.vcf.gz"
    bcftools sort "${DONOR_VCF_DIR}/tmp.biallelic_snps.vcf.gz" -Oz -o "$BIALLELIC_SNPS_VCF"
    tabix -p vcf "$BIALLELIC_SNPS_VCF"
  fi
```

### 4. Peak Restriction & Header Reordering

Only SNPs that overlap with accessible chromatin regions (ATAC peaks) are informative for this assay. We restrict the VCF to these coordinates. Finally, we deconstruct and rebuild the VCF header so the ##contig entries perfectly match the BAM chromosome order extracted in Step 0.

```R
# STEP 2: Restrict VCF to ATAC peak regions
  PEAKS_BED="${BED_FOLDER}/${SAMPLE_ID}/peaks.bed.gz"
  FILTERED_VCF="${DONOR_VCF_DIR}/donors.atac_filtered.biallelic_snps.vcf.gz"
  
  if [ ! -f "$FILTERED_VCF" ]; then
    log "Restricting VCF to ATAC peak regions..."
    bcftools view --threads "$THREADS" -R "$PEAKS_BED" "$BIALLELIC_SNPS_VCF" -Oz -o "${DONOR_VCF_DIR}/tmp.atac.vcf.gz"
    bcftools sort "${DONOR_VCF_DIR}/tmp.atac.vcf.gz" -Oz -o "$FILTERED_VCF"
    tabix -p vcf "$FILTERED_VCF"
  fi

  # STEP 3: Reorder VCF ##contig header entries
  REORDERED_VCF="${DONOR_VCF_DIR}/donors.atac_filtered.biallelic_snps.reordered.vcf.gz"
  if [ ! -f "$REORDERED_VCF" ]; then
    log "Reordering VCF ##contig header..."
    
    # Split header components
    bcftools view -h "$FILTERED_VCF" \vert{} grep '^##contig' > "${DONOR_VCF_DIR}/vcf.contigs.all.txt"
    bcftools view -h "$FILTERED_VCF" \vert{} grep '^##' \vert{} grep -v '^##contig' > "${DONOR_VCF_DIR}/vcf.header.meta.txt"
    bcftools view -h "$FILTERED_VCF" \vert{} grep '^#CHROM' > "${DONOR_VCF_DIR}/vcf.chromline.txt"

    # Re-build in BAM order
    : > "${DONOR_VCF_DIR}/vcf.contigs.ordered.txt"
    while read -r CHR; do
      grep -m1 -E "ID=${CHR}([,>])" "${DONOR_VCF_DIR}/vcf.contigs.all.txt" >> "${DONOR_VCF_DIR}/vcf.contigs.ordered.txt" || true
    done < "$BAM_ORDER_FILE"

    # Reassemble and zip
    bgzip -cd "$FILTERED_VCF" \vert{} awk '!/^#/' > "${DONOR_VCF_DIR}/vcf.body.tmp"
    cat "${DONOR_VCF_DIR}/vcf.header.meta.txt" "${DONOR_VCF_DIR}/vcf.contigs.ordered.txt" \
        "${DONOR_VCF_DIR}/vcf.chromline.txt" "${DONOR_VCF_DIR}/vcf.body.tmp" \vert{} bgzip -c > "$REORDERED_VCF"
    tabix -p vcf "$REORDERED_VCF"
  fi
```

### 5. Executing Demuxlet

With the VCF perfectly formatted and trimmed, demuxlet can rapidly assign cells to donors based on the CB (Cell Barcode) tags embedded in the BAM.

```R
# STEP 4: Assign cells to donors
  if [ ! -f "${DEMUXLET_DIR}/demuxlet.best" ]; then
    log "Running demuxlet..."
    popscle demuxlet \
      --tag-group CB \
      --sam "$BAM_FILE" \
      --vcf "$REORDERED_VCF" \
      --field GT \
      --sm-list "$DONOR_LIST_FILE" \
      --out "${DEMUXLET_DIR}/demuxlet"
    log "demuxlet finished for sample: $SAMPLE_ID"
  fi
done

log "=== scATAC demuxlet pipeline finished ==="
```
<FollowUp label="Want to parse the demuxlet results?" query="How should I parse the demuxlet.best output file to identify the singlet and doublet assignments for my downstream Seurat object?"/>