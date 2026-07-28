---
title: "JAGUAR-GRN"
date: 2026-07-26
draft: false
status: "in_progress"
summary: "A single-cell ATAC-seq pipeline designed to build a Gene Regulation Network to unveil the regulatory systems on JAGUAR Project Data."
tags: ["scATAC-seq", "JAGUAR", "R", "Python", "Seurat"]
github: "https://github.com/JAGUAR-LATAM/CCAN-Visualizer/"
stack: ["Snakemake", "R", "Seurat", "SLURM"]
---

## The Problem

Our existing pipeline establishes a strong foundation, but it leaves deeper biological insights hidden in the tables. By moving toward Gene Regulation Networks (GRNs), we can understand the *behavior* of the data.

| Existing Capabilities | New GRN Capabilities |
|---|---|
| Filter scATAC-seq data by quality | Discover synchronously activated regions |
| Group cells by accessible regions | Map enhancer-to-target interactions |
| Perform differential accessibility | Identify regulatory topics and cellular archetypes |

My work in this phase is to unveil how these patterns manifest and verify whether they align with known biological regulatory systems.

---

## Approach: Cicero & Co-accessibility

To build the main highway of our network, I use **Cicero**. Cicero reads scATAC-seq data to identify when peaks are co-accessible (accessible simultaneously). Because our data lacks a time-series dimension, we calculate a consensus across all cells.

This allows us to construct **CCANs** (Cis-Co-Accessibility Networks) — clusters of near-proximity peaks that form interacting blocks, rather than just simple pairwise interactions. 

### Pipeline Flow
`[Input: Seurat .rds]` -> 

`[Monocle3 CDS]` -> 

`[Cicero Co-accessibility]` -> 

`[CCAN Modules]` -> 

`[Gene Activity Scores]` -> 

`[Output: Integrated Seurat]`

---

## Implementation Guide

Below is the step-by-step R script to execute this pipeline. It handles argument parsing for cluster environments (like SLURM|LSF), calculates co-accessibility, and integrates the results back into a Seurat object.

### 1. Environment Setup and Argument Parsing
First, we load the required libraries and parse command-line arguments to avoid hardcoding file paths.

Here is the usage in case you're code-blinded:

```bash
Rscript script_name.R <path_to_seurat.rds> <output_prefix> [project_root_optional]
# Example
Rscript create_cicero_cds.r integrated_seurat.rds cicero_test ~/JAGUAR/GRN/
```

```R
options(repos = c(CRAN = "[https://cloud.r-project.org](https://cloud.r-project.org)"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script_name.R <path_to_seurat.rds> <output_prefix> [project_root_optional]", call. = FALSE)
}

seurat_object_path <- args[1]
prefix <- args[2]
project_root <- ifelse(length(args) >= 3, args[3], file.path(Sys.getenv("HOME"), "JAGUAR/GRN"))

# Package Loading (assumes dependencies are installed)
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SeuratWrappers)
  library(monocle3)
  library(cicero) 
  library(EnsDb.Hsapiens.v86)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

# Define Paths
output_dir <- file.path(project_root, paste0("3_output/cicero_", prefix))
figure_dir <- file.path(project_root, "4_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
```

### 2. Creating the Cicero CDS Object

Cicero relies on Monocle3's Cell Data Set (CDS) format. We check if an object already exists to save computation time; otherwise, we generate it from the input Seurat object.

```R
cicero_cds_path <- file.path(project_root, paste0("3_output/cicero_", prefix, "_cds"))

if (dir.exists(cicero_cds_path)) {
  print("[DEBUG] Loading existing Cicero object from disk...")
  cicero_cds <- load_monocle_objects(directory_path = cicero_cds_path)
  input_seurat <- readRDS(seurat_object_path)
  
} else {
  print("[DEBUG] Creating new Cicero object...")
  input_seurat <- readRDS(seurat_object_path)
  
  # Convert Seurat to Monocle3 CDS
  input_cds <- as.cell_data_set(x = input_seurat)
  input_cds <- cluster_cells(input_cds, reduction_method = "UMAP")
  
  # Extract UMAP and create Cicero CDS
  umap_coords <- reducedDims(input_cds)$UMAP
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords, k = 50)
  
  save_monocle_objects(cicero_cds, directory_path = cicero_cds_path)
}
```

### 3. Calculating Co-accessibility

Using the hg38 human genome reference, we calculate the connections between accessible peaks.

```R
genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
genome_df <- data.frame("chr" = names(genome), "length" = genome)

print("[DEBUG] Calculating co-accessibility connections...")
conns <- run_cicero(cicero_cds, genome_df, sample_num = 100)

# Export for external tools (e.g., CellOracle)
conns_filtered <- subset(conns, coaccess > 0)
coaccess_csv_path <- file.path(output_dir, paste0(prefix, "_coaccess_cicero.csv"))
write.csv(conns_filtered, coaccess_csv_path, row.names = FALSE, quote = FALSE)
```

### 4. Annotation and Visualization

We fetch human gene annotations to map our peaks to known genes, then generate a visualization of a specific locus (e.g., OAS1).


```R
print("[DEBUG] Fetching human gene annotations...")
tx <- transcripts(EnsDb.Hsapiens.v86, columns = c("gene_id", "tx_id", "gene_name"))
gene_anno <- as.data.frame(tx)

# Match column names for Cicero
gene_anno$chromosome <- paste0("chr", gene_anno$seqnames)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$tx_id
gene_anno$symbol <- gene_anno$gene_name

pdf_filename <- paste0(prefix, "_Cicero_Connections_OAS1_Locus.pdf")
pdf(file.path(figure_dir, pdf_filename), width = 40, height = 20)
plot_connections(conns, "chr12", 112490000, 113200000, # <- Around OAS1 locus
                 gene_model = gene_anno,
                 coaccess_cutoff = 0.1,
                 connection_width = 0.5,
                 alpha_by_coaccess = TRUE,
                 collapseTranscripts = "longest")
dev.off()
```

### 5. Identifying CCANs

We group the co-accessible links into Louvain communities (CCANs) to find broader regulatory blocks.

```R
print("[DEBUG] Searching for CCAN modules...")
CCAN_assigns <- generate_ccans(conns_filtered)

ccans_csv_path <- file.path(output_dir, paste0(prefix, "_CCANs.csv"))
write.csv(CCAN_assigns, ccans_csv_path, row.names = FALSE)
```

### 6. Calculating Gene Activity Scores

By linking distal regulatory elements to local promoters, we estimate the resulting transcriptional activity (translating ATAC signals to expected RNA expression).

```R
print("[DEBUG] Calculating Gene Activity Scores...")

# Extract Transcription Start Sites (TSS)
pos <- subset(gene_anno, strand == "+")
pos <- pos[order(pos$start), ] 
pos <- pos[!duplicated(pos$transcript), ] 
pos$end <- pos$start + 1 

neg <- subset(gene_anno, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE), ] 
neg <- neg[!duplicated(neg$transcript), ] 
neg$start <- neg$end - 1

gene_annotation_sub <- rbind(pos, neg)[, c("chromosome", "start", "end", "symbol")]
names(gene_annotation_sub)[4] <- "gene"

# Annotate CDS and build unnormalized matrix
if (!exists("input_cds")) input_cds <- as.cell_data_set(x = input_seurat)
input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)
unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, !Matrix::colSums(unnorm_ga) == 0]

# Normalize Matrix
if (is.null(pData(input_cds)$num_genes_expressed)) {
  num_genes <- ifelse(!is.null(pData(input_cds)$nFeature_ATAC), 
                      pData(input_cds)$nFeature_ATAC, 
                      Matrix::colSums(counts(input_cds) > 0))
} else {
  num_genes <- pData(input_cds)$num_genes_expressed
}
names(num_genes) <- row.names(pData(input_cds))

cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)

matrix_rds_path <- file.path(output_dir, paste0(prefix, "_gene_activities_matrix.rds"))
saveRDS(cicero_gene_activities, matrix_rds_path)
```

### 7. Seurat Integration

Finally, we load the Gene Activity matrix back into the original Seurat object as a new assay, allowing us to map activity onto our standard UMAPs.

```R
print("[DEBUG] Integrating Gene Activity Scores into Seurat...")
cicero_gene_activities <- readRDS(matrix_rds_path)
cicero_gene_activities <- cicero_gene_activities[, colnames(input_seurat)]

# Create and set the new Assay
input_seurat[["ACTIVITY"]] <- CreateAssayObject(counts = cicero_gene_activities)
DefaultAssay(input_seurat) <- "ACTIVITY"

# Visualize
FeaturePlot(input_seurat, features = c("MNT", "TTC7A", "CARS2"), pt.size = 0.5)

# Verify integration
matches <- grep("CDC", rownames(cicero_gene_activities), ignore.case = TRUE, value = TRUE) 
print(matches)
```

## NEXT-STEPS: pycisTopic

To improve the confidence in our results, we can validate our findings using **pycisTopic**. This tool uses Latent Dirichlet Allocation (LDA) to group co-accessible chromatin regions into distinct regulatory topics and classify cell states.

### 1. Data Preparation: "Recoger nuestras chivas"

Before diving into the Python environment, we need to get our files in order. `pycisTopic` relies on a core data structure: a dictionary that links sample IDs to their corresponding ATAC fragment files. 

Because our data currently lives inside a Seurat object, we need to extract four specific deliverables to feed into the Python pipeline:

| Output File | Purpose |
|---|---|
| `integrated_peaks.bed` | The genomic coordinates of all accessible regions. |
| `valid_barcodes.txt` | The filtered list of cells that passed Seurat's quality control. |
| `sample_fragment_mapping.tsv` | The key-value pairing of datasets to their cluster fragment paths. |
| `seurat_metadata_with_fragments.tsv` | The full metadata table, updated to include the exact fragment path for each cell. |

---

### Implementation Guide

The following R script extracts these required components directly from the integrated Seurat object.

#### Step 1: Environment Setup
First, we load the object and define our output directories on the cluster.

> PERSONAL NOTE: STOP HARDCODING VARIABLES USE ARGs INSTEAD!!!!!!

```R
library(Seurat)
library(Signac)

# Load Seurat object
seurat_obj <- readRDS("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/de8/yascp_scatac/Downsampling_analysis/complete_analysis/outs/integrated_filtered.rds")

outdir <- "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jv8/GRN/2_data"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Loading Seurat object... OK\n")
cat("Output directory:", outdir, "\n\n")
```

### Step 2: Export Peaks and Valid Barcodes

We extract the GRanges object containing our peaks and convert it to a standard BED format, followed by exporting the names of the cells that passed filtering.

```R
cat("=== EXPORTING PEAKS ===\n")
peaks_gr <- granges(seurat_obj[["ATAC"]])
peaks_df <- as.data.frame(peaks_gr)[, c("seqnames", "start", "end")]
write.table(peaks_df, file.path(outdir, "integrated_peaks.bed"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("OK: Peaks exported\n\n")

cat("=== EXPORTING BARCODES ===\n")
valid_barcodes <- colnames(seurat_obj)
writeLines(valid_barcodes, file.path(outdir, "valid_barcodes.txt"))
cat("OK: Exported", length(valid_barcodes), "barcodes\n\n")
```

### Step 3: Extract and Verify Fragment Paths

We pull the fragment paths directly from the Seurat ATAC assay and verify that the files actually exist on the lustre file system.

```R
cat("=== EXTRACTING FRAGMENT PATHS ===\n")
frag_list <- Fragments(seurat_obj[["ATAC"]])
cat("Number of fragment files:", length(frag_list), "\n")

# Extract paths
frag_paths <- sapply(frag_list, function(f) f@path)
names(frag_paths) <- NULL

cat("=== VERIFYING FILES EXIST ===\n")
for (i in seq_along(frag_paths)) {
  exists <- file.exists(frag_paths[i])
  status <- ifelse(exists, "✓", "✗")
  cat(status, "Fragment", i, ":", frag_paths[i], "\n")
}
cat("\n")
```

### Step 4: Create the Sample Mapping Dictionary

We map the unique datasets found in the metadata to their corresponding fragment files. This will be the foundation of the pycisTopic dictionary.

```R
cat("=== CREATING SAMPLE MAPPING ===\n")

unique_datasets <- sort(unique(seurat_obj$dataset[!is.na(seurat_obj$dataset)]))

if (length(unique_datasets) != length(frag_paths)) {
  cat("WARNING: Number of datasets (", length(unique_datasets), 
      ") does not match number of fragment files (", length(frag_paths), ")\n")
}

mapping_df <- data.frame(
  sample_id = seq_along(frag_paths),
  dataset = unique_datasets[seq_along(frag_paths)],
  fragment_path = frag_paths,
  stringsAsFactors = FALSE
)

write.table(mapping_df, file.path(outdir, "sample_fragment_mapping.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("OK: Mapping saved to sample_fragment_mapping.tsv\n\n")
```

### Step 5: Update and Export Metadata

Finally, we append the specific fragment path to every individual cell in the metadata and export the complete table.

```R
cat("=== ADDING SAMPLE INFO TO METADATA ===\n")

# Map dataset to fragment_path using the mapping_df
seurat_obj$fragment_path <- mapping_df$fragment_path[match(seurat_obj$dataset, mapping_df$dataset)]

cat("Fragment path assignment check:\n")
print(table(!is.na(seurat_obj$fragment_path)))
cat("\n")

cat("=== EXPORTING METADATA ===\n")
write.table(seurat_obj@meta.data, file.path(outdir, "seurat_metadata_with_fragments.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
cat("OK: Full metadata exported\n")
```

## FUTURE WORK

Now I will try to complete pycistopic tutorial. At this point if you have a seurat object with the previous data you could follow it too. btw cloud computing sucks buy me a home lab plspls #grant #jovenpromesa