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

That little diagram above is really just step 1. In practice all of the Cicero-branch scripts get chained together by one driver, `1_scripts/cicero/run_pipeline.sh`, which submits each step to LSF and wires a `-w "done(...)"` dependency onto whatever it needs finished first — so the whole 6-step branch can be launched with a single `./run_pipeline.sh` and left to run unattended (steps 3-6 for one prefix take minutes-to-hours; step 1 is the multi-day one). Here's the full chain, with a checkbox for whether I've actually written up that step in this post yet:

- [x] **Step 1 - `create_cicero_cds.r`**: Seurat → Cicero CDS → co-accessibility → CCANs → gene-activity matrix. The long one (~48h budget, only resumable up to the cached CDS). Walked through below in ["Create a cds object and Run Cicero"](#create-a-cds-object-and-run-cicero).
- [x] **Step 2 - `TestingCCANSandCONNS.R`**: generates the per-CCAN track-plot PDFs (peaks + arcs + genes). Runs in parallel with steps 3-6, gated only on step 1. Walked through below in ["Analyze and Plot each CCANs"](#analyze-and-plot-each-ccans).
- [ ] **Step 3 - `assign_ccan_genes.r`**: assigns genes to each CCAN via GREAT's regulatory-domain lookup (cheap, no statistical test). *Working on it...*
- [ ] **Step 4 - `annotate_ccans_go.r`** (LSF array, one task per CCAN stripe) + merge: runs GREAT GO-term enrichment per CCAN (GO:BP by default, ~80-100s/CCAN — the expensive one) and merges the array's part-files back into one CSV. *Not written up yet.*
- [ ] **Step 5 - `cluster_ccans_by_go.r`**: builds a CCAN × GO-term TF-IDF matrix from step 3's output and runs PCA → UMAP → Leiden clustering, validated against genomic locality and GO/immune-keyword ground truth. Needs steps 3 and 4. *Not written up yet.*
- [ ] **Step 6 - `build_ccan_index_go.r`**: assembles the final per-CCAN JSON manifest (region/size/PDF/GO/cluster fields combined) for an external visualizer. Needs steps 3 and 5. *Not written up yet.*
- [ ] **Optional - RegulomeDB branch** (`annotate_ccans_regulomedb.r` LSF array + `merge_regulomedb_shards.r`): off by default (`--run-regulomedb`); queries the RegulomeDB REST API per peak for independent chromatin/TF-binding evidence. Runs in parallel with steps 3-6, gated only on step 1. Not yet wired into step 6's JSON either way — still future work. *Not written up yet.*

So this post currently covers steps 1 and 2 of the six — steps 3-6 (plus the optional RegulomeDB branch) are next on the list to document.

---

## Create a cds object and Run Cicero

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
  
}
```

If it didn't exist, then we we need to create it using a cell_data_set, but I don't have that, I only have a seurat object, so using a function from SeuratWrapper called `as.cell_data_cell()` (To see more about cell_data_set class please see: https://cole-trapnell-lab.github.io/monocle3/docs/getting_started/#cell_data_set). 

```R
 else {
  print("[DEBUG] Creating new Cicero object...")
  input_seurat <- readRDS(seurat_object_path)
  
  # Convert Seurat to Monocle3 CDS
  input_cds <- as.cell_data_set(x = input_seurat)
```

Then to run cicero with the cds object first I need to cluster the cells using a reduction method in this case UMAP, in other words we use all the data from the cell accessibility and give them cordinates given this accessibility data, as If the data give a cell a place to stay. then extract the cordinates from the cells and give that to the function `make_cicero_cds()` and save the object using `save_monocle_objects()`. Note: by default make_cicero_cds function k is equal to 50, k is the number of nearest neighbors used to build overlapping metacells.

```R
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
```

Then I filter the connections to have a coaccessibility greater than 0, meaning that I only kept peaks that have a connection and save it into a csv (a table) to use it in other tools like `CellOracle`.

```R
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

By linking distal regulatory elements to local promoters, we estimate the resulting transcriptional activity (translating ATAC signals to expected RNA expression). First I need to get the gene that is close to the peak, that is Gene Activity, like a pseudo RNA-seq but instead of having transcriptional data we have the accessibility of the DNA, so we guess the transcriptional activity.

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

## Analyze and Plot each CCANs

### 1. Load data

Now that we have the results of running cicero it is time to see what is inside.
First we load our packages.

```r
# Script: TestingCCANSandCONNS.R
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SeuratWrappers)
  library(monocle3)
  library(cicero)
  library(EnsDb.Hsapiens.v86)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(VariantAnnotation)
  library(igraph)
  library(dplyr)
  library(ggforce)
  library(ggplot2)
  library(patchwork)  # For combining plots (plot_layout(), `/` operator)
  library(scales)     # For rescale()
})
```

The argumentss mirrors `create_cicero_cds.r's <output_prefix> [project_root] `convention:
that script names both the output directory and every file inside it after
the same `prefix` (3_output/cicero_<prefix>/<prefix>_CCANs.csv, etc).
``` bash
#Usage: 
Rscript TestingCCANSandCONNS.R [dir_prefix] [file_prefix] [project_root]
```

```r
args <- commandArgs(trailingOnly = TRUE)

dir_prefix   <- if (length(args) >= 1) args[1] else "full_cicero"
file_prefix  <- if (length(args) >= 2) args[2] else "cells"
project_root <- if (length(args) >= 3) args[3] else file.path(Sys.getenv("HOME"), "JAGUAR/GRN")

output_dir <- file.path(project_root, paste0("3_output/cicero_", dir_prefix))
figure_dir <- file.path(project_root, "4_figures")

conns <- read.csv(file.path(output_dir, paste0(file_prefix, "_coaccess_cicero.csv")))
ccans <- read.csv(file.path(output_dir, paste0(file_prefix, "_CCANs.csv")))
```


### 2. Descriptive statistics


Then we check our numbers to validate our data:

- How Many CCANs Exists?
- How Many Peaks per CCANs Exists?

We should expect that a CCAN have at least 2 peaks.

```r

print(head(ccans))

cat("||||How Many CCANs Exists?||||")
cat("There are",length(unique(ccans$CCAN)),"CCANs")

cat("||||How Many Peaks per CCANs Exists?||||")
hist(ccans$CCAN, breaks = length(unique(ccans$CCAN)))


empty <- c()

ccans_ids <- unique(ccans$CCAN)

region_df <- data.frame()
```

This Step creates a data frame per CCAN where only kept the start and end of the ccan instead of having all the peaks.

In english this chunk of code would be, for each CCAN take that ccan from the general dataset, then extract the position of the range of each peak from that CCAN
then create a data frame where you include the name of that CCAN, which chromosome belongs to, the minimum value of all the positions (~start~) and the maximum value of all the positions (~end~) and the number of peaks. 
```r


for(i in ccans_ids){

  sub <- subset(ccans, CCAN == i)

  coords <- do.call(
    rbind,
    strsplit(sub$Peak, "-")
  )

  region_df <- rbind(
    region_df,
    data.frame(
      ccan = i,
      chromosome = unique(coords[,1]),
      start = min(as.numeric(coords[,2])),
      end = max(as.numeric(coords[,3])),
      n_peaks = nrow(sub)
    )
  )
}
```

With the new data frame called `region_df`, we calculate the size of each CCAN by substracting the end and the start (the distance lmao).

Then, make a plot comparing the size of each CCAN and the number of peaks of each CCAN to answer **Do bigger CCANs contain more peaks?**, The plot confirmed it (I don't have the plot lmao but trust me (that isn't a good scientific practice but do it))

Also I checked the top 20 ccans by number of peaks and by size to see if they're the same. (They were the same)
```r
region_df$size_bp <- region_df$end - region_df$start
region_df$size_kb <- region_df$size_bp / 1000
region_df$size_mb <- region_df$size_bp / 1e6

# Do bigger CCANs contain more peaks?

plot(
  region_df$size_kb,
  region_df$n_peaks,
  xlab = "Size (kb)",
  ylab = "N Peaks"
)

# Answer: Yes, Like it should?

# Top 20 ccans by number of peaks (picudos)
head(
  region_df[order(-region_df$n_peaks), ],
  20
)

# Top 20 ccans by size (grandotes)
head(
  region_df[order(-region_df$size_bp), ],
  20
)
```

### 3. Annotate CCANs

Now it's time to see which genes are nerby the peaks of the CCANs, We get the information of each gene by giving the start and the end to `GRanges`.

```r

gr_ccan <- GRanges(
  seqnames = region_df$chromosome,
  ranges = IRanges(
    start = region_df$start,
    end = region_df$end
  )
)

gene_gr <- genes(EnsDb.Hsapiens.v86)
```
We got the gene annotation from `EnsDb.Hsapiens.v86`, to then find the overlaps using the function `FindOverlaps()` between our ranges and the genes from the database.

> `gene_gr` is reused below both by the single-CCAN worked example (Section 4) and by the generalized plot_ccan_track() function (Section 5) -- computed once here rather than re-derived in each place.
```r
seqlevelsStyle(gene_gr) <- "UCSC"

hits <- findOverlaps(
  gr_ccan,
  gene_gr
)

ccan_genes <- data.frame(
  ccan = region_df$ccan[queryHits(hits)],
  gene = gene_gr$gene_name[subjectHits(hits)]
)

# Test to find OAS family genes

subset(
  ccan_genes,
  gene %in% c("OAS1","OAS2","OAS3","OASL")
)

# see how many genes are per ccan

genes_per_ccan <- aggregate(
  gene ~ ccan,
  ccan_genes,
  function(x) length(unique(x))
)

head(
  genes_per_ccan[
    order(-genes_per_ccan$gene),
  ],
  20
)
```
### 4.- An example using OAS 

Using the CCAN that contain OAS (`ccan_id = 1025`) we would process that info to plot it as a coaccessibility per peak of that CCAN. First make a graph of each conection from the CCAN.
```r

ccan_id <- 1025

# Subset the ccan that contain oas

oas_peaks <- subset(
  ccans,
  CCAN == ccan_id
)$Peak

oas_conns <- subset(
  conns,
  Peak1 %in% oas_peaks &
    Peak2 %in% oas_peaks
)

# Create a graph of the connexions

g <- graph_from_data_frame(
  oas_conns[, c("Peak1","Peak2")],
  directed = FALSE
)

deg <- degree(g)

hubs <- names(sort(deg, decreasing = TRUE)[1:20])
```

Now that we have our top 20 "hub" peaks (the ones with the most connections in the graph), we turn them into an actual data frame so we can eventually plot them: split each peak ID (`chr-start-end`) apart and compute its midpoint (`center`) — that midpoint is the x-coordinate we'll later hand to ggplot to draw something at that spot.

Then, separately, we score peaks a different way: instead of just counting connections (degree), we add up how *strong* those connections are. For every peak, sum the `coaccess` value of every connection it's part of — that gives `total_strength`. Taking the top 10 peaks by that score gives us a coaccess-weighted hub list, which is a bit more informative than raw degree alone.

```r
# Make a data frame that contain the center of each peak

hub_df <- do.call(
  rbind,
  lapply(hubs, function(x){

    y <- strsplit(x, "-")[[1]]

    data.frame(
      chr = y[1],
      start = as.numeric(y[2]),
      end = as.numeric(y[3]),
      center = (as.numeric(y[2]) + as.numeric(y[3]))/2
    )
  })
)

# Given the overall coaccess per peak sum all that and obtain the total strength of the node

node_strength <- bind_rows(
  oas_conns %>%
    select(Peak = Peak1, coaccess),

  oas_conns %>%
    select(Peak = Peak2, coaccess)
) %>%
  group_by(Peak) %>%
  summarise(
    degree = n(),
    mean_coaccess = mean(coaccess),
    total_strength = sum(coaccess)
  )

node_strength |> arrange(desc(total_strength))

top10 <- node_strength |>
  arrange(desc(total_strength)) |>
  head(10)
```

Quick gut-check before going further: I hardcoded the coordinates of the single biggest hub peak I'd already spotted (`chr12:113040427-113041384`) and just asked, in the console, "what genes are actually near this thing?" using `findOverlaps()` and `distanceToNearest()`. Nothing here gets saved — it's scratch code, not part of the real pipeline — but it's how I first sanity-checked that these hub peaks were landing near genes I recognized before trusting the rest of the analysis.

```r
peak_gr_peak_top <- GRanges(
  "chr12",
  IRanges(
    113040427,
    113041384
  )
)

findOverlaps(
  peak_gr_peak_top,
  gene_gr
)

dists <- distanceToNearest(
  peak_gr_peak_top,
  gene_gr,
  select = "all"
)
```

Now for real, across all 10 hub peaks at once: turn them into a proper `GRanges` object, then ask two questions — what's the single nearest gene to each hub in general, and specifically, how close is each hub to the known OAS genes (`OAS1`, `OAS2`, `OAS3`)? That second question is the whole point of this worked example — we already suspect this CCAN regulates the OAS locus, so we want to see the hubs "reaching toward" it.

```r
peaks_split <- do.call(rbind, strsplit(top10$Peak, "-"))

hubs_gr <- GRanges(
  seqnames = peaks_split[,1],
  ranges = IRanges(start = as.numeric(peaks_split[,2]),
                   end = as.numeric(peaks_split[,3]))
)

top10_nearest_hits <- distanceToNearest(
  hubs_gr,
  gene_gr
)

oas_hits <- distanceToNearest(
  hubs_gr,
  subset(
    gene_gr,
    gene_name %in% c("OAS1","OAS2","OAS3")
  ),
  select = "all"
)

oas_conns |>
  arrange(desc(coaccess)) |>
  head(20)
```

(That last one's just me eyeballing the 20 strongest connections in the console — no assignment, no save, just looking.)

Time to put all of this together into one clean table. We pull the gene name and distance out of the `distanceToNearest()` result from before, and build `summary_table`: one row per hub peak, with its connection degree, mean/total coaccess strength, and its nearest gene + how far away that gene is.

```r
nearest_gene <- gene_gr$gene_name[
  subjectHits(top10_nearest_hits)
]

nearest_distance <- mcols(top10_nearest_hits)$distance

summary_table <- data.frame(
  Peak = top10$Peak,
  Degree = top10$degree,
  Mean_Coaccess = round(top10$mean_coaccess,3),
  Total_Strength = round(top10$total_strength,3),
  Nearest_Gene = nearest_gene,
  Distance = nearest_distance
)

summary_table
```

Not done yet — let's enrich `summary_table` with a few OAS-specific columns: does each hub peak physically overlap an OAS gene body, how far is it from the OAS locus as a whole (treating `OAS1`/`OAS2`/`OAS3` as one combined region), and a simple `HubScore` (degree × mean coaccess) to rank hubs by "how connected and how strong" in one number.

```r
coords <- do.call(
  rbind,
  strsplit(summary_table$Peak, "-")
)

summary_table$chr <- coords[,1]
summary_table$start <- as.numeric(coords[,2])
summary_table$end <- as.numeric(coords[,3])

oas_gr <- subset(
  gene_gr,
  gene_name %in% c("OAS1","OAS2","OAS3")
)

hits_oas <- overlapsAny(
  hubs_gr,
  oas_gr
)

summary_table$Overlaps_OAS <- hits_oas

oas_locus <- GRanges(
  "chr12",
  IRanges(
    start = min(start(oas_gr)),
    end = max(end(oas_gr))
  )
)

dist_oas <- distanceToNearest(
  hubs_gr,
  oas_locus
)

summary_table$Distance_to_OAS <- mcols(dist_oas)$distance

summary_table$HubScore <-
  summary_table$Degree *
  summary_table$Mean_Coaccess
```

With the exploration done, it's time to actually build the plot. First step: don't try to draw *every* connection in the CCAN, that's way too noisy — only keep the strongest ones. We take everything at or above the 95th percentile of coaccess (`threshold`), drop duplicate pairs (a connection A-B is the same as B-A, so we only keep rows where `Peak1 < Peak2` alphabetically), and compute each connection's midpoint-to-midpoint x-coordinates (`x1`, `x2`) — that's what lets us draw an arc from one peak to the other. `curvature`, `alpha_val`, and `width_val` are just the connection's coaccess strength rescaled into plotting ranges, so stronger connections get taller, more opaque, thicker arcs.

```r
threshold <- quantile(
  oas_conns$coaccess,
  0.95
)

oas_sig <- subset(
  oas_conns,
  coaccess >= threshold
)

oas_sig <- oas_sig[oas_sig$Peak1 < oas_sig$Peak2, ]

peak_center <- function(peak) {

  x <- strsplit(peak, "-")[[1]]

  start <- as.numeric(x[2])
  end   <- as.numeric(x[3])

  (start + end) / 2
}

arc_df <- oas_sig

arc_df$x1 <- sapply(
  arc_df$Peak1,
  peak_center
)

arc_df$x2 <- sapply(
  arc_df$Peak2,
  peak_center
)

arc_df$curvature <- rescale(
  arc_df$coaccess,
  to = c(0.1, 0.8)
)

arc_df$alpha_val <- rescale(arc_df$coaccess, to = c(0.45, 0.95))
arc_df$width_val <- rescale(arc_df$coaccess, to = c(0.8, 2.8))
```

Now the actual plot. The idea is a genome-browser-style track: three panels stacked on top of each other, all sharing the exact same x-axis (genomic position), so everything lines up vertically. From top to bottom: the co-accessibility arcs, the peaks themselves, and the genes underneath — so you can visually trace "this arc connects this peak to that peak, and oh look, that peak sits right on top of this gene." We build each panel as its own ggplot object and glue them together with `patchwork` (the `/` operator) at the end.

```r
# --------------------------------------------------------------------------
# Plot coaccessibility + peaks + genes
# --------------------------------------------------------------------------
# Panel 1 (top)  -> Arch
# Panel 2 (mid)   -> Peaks
# Panel 3 (below)   -> Genes

# 0. Whole Window Size

ccan_region <- subset(region_df, ccan == ccan_id)

region_chr   <- ccan_region$chromosome
region_start <- ccan_region$start - 5000
region_end   <- ccan_region$end + 5000

# 1. Peak Panel Plot

peaks_split_all <- do.call(rbind, strsplit(oas_peaks, "-"))

peaks_df <- data.frame(
  chr   = peaks_split_all[, 1],
  start = as.numeric(peaks_split_all[, 2]),
  end   = as.numeric(peaks_split_all[, 3])
)

peaks_df$width  <- peaks_df$end - peaks_df$start
peaks_df$is_hub <- oas_peaks %in% top10$Peak

peak_plot <- ggplot(peaks_df) +
  geom_rect(
    aes(xmin = start, xmax = end, ymin = 0.2, ymax = 0.8, fill = is_hub)
  ) +
  scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), guide = "none") +
  coord_cartesian(xlim = c(region_start, region_end)) +
  labs(title = "ATAC peaks (red = top10 hubs)") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 2)))

# 2. Arch Plot Panel

arc_plot <- ggplot() +
  coord_cartesian(xlim = c(region_start, region_end), ylim = c(0, 1))

for (i in seq_len(nrow(arc_df))) {

  arc_plot <- arc_plot +
    geom_curve(
      data = arc_df[i, ],
      aes(x = x1, y = 0, xend = x2, yend = 0),
      curvature  = -arc_df$curvature[i],
      alpha      = arc_df$alpha_val[i],
      linewidth  = arc_df$width_val[i],
      color      = "steelblue4"
    )
}

arc_plot <- arc_plot +
  labs(title = paste0("Co-Accsessibility (>= cuantil 0.95, coaccess >= ",
                      round(threshold, 3), ")")) +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 2)))

# 3. Gene Plot Panel

genes_in_region <- subset(
  gene_gr,
  as.character(seqnames(gene_gr)) == region_chr &
    start(gene_gr) < region_end &
    end(gene_gr)   > region_start
)

genes_df <- as.data.frame(genes_in_region)
genes_df$is_oas <- genes_df$gene_name %in% c("OAS1", "OAS2", "OAS3")

genes_df$gene_name <- factor(
  genes_df$gene_name,
  levels = unique(
    genes_df$gene_name[
      order(genes_df$start)
    ]
  )
)

genes_plus  <- subset(genes_df, strand == "+")
genes_minus <- subset(genes_df, strand == "-")

gene_plot <- ggplot() +
  { if (nrow(genes_plus) > 0)
    geom_segment(
      data = genes_plus,
      aes(x = start, xend = end, y = gene_name, yend = gene_name, color = is_oas),
      linewidth = 3,
      arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "open")
    )
  } +
  { if (nrow(genes_minus) > 0)
    geom_segment(
      data = genes_minus,
      aes(x = start, xend = end, y = gene_name, yend = gene_name, color = is_oas),
      linewidth = 3,
      arrow = arrow(length = unit(0.2, "cm"), ends = "first", type = "closed")
    )
  } +
  scale_color_manual(values = c(`FALSE` = "grey40", `TRUE` = "firebrick"), guide = "none") +
  coord_cartesian(xlim = c(region_start, region_end)) +
  labs(
    title = "Genes in Region",
    x = paste0("Position: ", region_chr, " (bp)")
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 2)),
    axis.title.y = element_blank()
  )

# 4. Combine Plots

final_plot <- arc_plot / peak_plot / gene_plot +
  plot_layout(heights = c(2, 0.6, 1.8))

final_plot

ggsave(
  filename = file.path(figure_dir, paste0("CCAN_", ccan_id, "_track_plot.png")),
  plot = final_plot,
  width = 9, height = 6.5, dpi = 300
)
```

### 5. Generalizing the Plot into a Reusable Function

Everything above was hardcoded for one CCAN (`ccan_id <- 1025`, the OAS one). Obviously I don't want to copy-paste that whole block by hand for every single CCAN in the dataset — so this section just wraps the exact same recipe (find hubs, filter by significance, build the 3-panel plot, save it) into one function, `plot_ccan_track()`, that takes a `ccan_id` and does the rest automatically. It's basically Section 4 again, but parameterized and with some extra guard rails added for the cases that didn't come up with our one hand-picked OAS example:

- If a CCAN has zero connections at all, don't try to plot it — just skip it and say so.
- If the 95th-percentile threshold happens to filter out *everything* (can happen with very few connections), fall back to just showing the top 10 strongest ones instead of crashing.
- If there's only one significant connection, `rescale()` has nothing to rescale against, so we hardcode sensible default curvature/alpha/width values instead.
- Output files are saved into a per-chromosome, per-`dir_prefix` folder (e.g. `4_figures/full_cicero/chr12/CCAN_1025_track_plot.pdf`) instead of one flat folder — see the comment in the code below for why: CCAN IDs are just local cluster labels, not stable identifiers, so "CCAN 42" from one run and "CCAN 42" from another run (say, a different cell type) are unrelated regions that would otherwise silently overwrite each other's plot.
- The gene panel's height is now dynamic (`gene_panel_weight`) instead of fixed, so a CCAN overlapping 30 genes doesn't squash them all unreadably into the same space as a CCAN overlapping 2.

```r
# Define the Generalized Plotting Function
# (gene_gr was already prepared in Section 3 and is passed in as an argument)
plot_ccan_track <- function(ccan_id, ccans_df, conns_df, region_df, gene_gr, fig_dir) {

  # --- Subsetting Data for Current CCAN ---
  current_peaks <- subset(ccans_df, CCAN == ccan_id)$Peak
  current_conns <- subset(conns_df, Peak1 %in% current_peaks & Peak2 %in% current_peaks)

  # Safety check: If no connections exist, skip plotting
  if(nrow(current_conns) < 1) {
    message(paste("Skipping CCAN", ccan_id, "- No connections found."))
    return(NULL)
  }

  # --- Calculate Hubs ---
  # Total strength of the node
  node_strength <- bind_rows(
    current_conns %>% select(Peak = Peak1, coaccess),
    current_conns %>% select(Peak = Peak2, coaccess)
  ) %>%
    group_by(Peak) %>%
    summarise(
      degree = n(),
      mean_coaccess = mean(coaccess),
      total_strength = sum(coaccess)
    )

  top10 <- node_strength %>% arrange(desc(total_strength)) %>% head(10)

  # Extract hubs GRanges to find nearest genes (to highlight them dynamically)
  peaks_split <- do.call(rbind, strsplit(top10$Peak, "-"))
  hubs_gr <- GRanges(
    seqnames = peaks_split[,1],
    ranges = IRanges(start = as.numeric(peaks_split[,2]),
                     end = as.numeric(peaks_split[,3]))
  )

  # Find genes nearest to the top 10 hubs to highlight them in the plot
  nearest_hits <- distanceToNearest(hubs_gr, gene_gr)
  hub_genes_to_highlight <- unique(gene_gr$gene_name[subjectHits(nearest_hits)])

  # --- Filter Connections by Significance ---
  threshold <- quantile(current_conns$coaccess, 0.95, na.rm = TRUE)
  current_sig <- subset(current_conns, coaccess >= threshold)

  # Safety check: If threshold is too high and removes everything, fallback to top connections
  if(nrow(current_sig) == 0) current_sig <- head(arrange(current_conns, desc(coaccess)), 10)

  current_sig <- current_sig[current_sig$Peak1 < current_sig$Peak2, ]

  peak_center <- function(peak) {
    x <- strsplit(peak, "-")[[1]]
    (as.numeric(x[2]) + as.numeric(x[3])) / 2
  }

  arc_df <- current_sig
  arc_df$x1 <- sapply(arc_df$Peak1, peak_center)
  arc_df$x2 <- sapply(arc_df$Peak2, peak_center)

  # Handle scaling safely if there's only 1 row
  if(nrow(arc_df) > 1) {
    arc_df$curvature <- rescale(arc_df$coaccess, to = c(0.1, 0.8))
    arc_df$alpha_val <- rescale(arc_df$coaccess, to = c(0.45, 0.95))
    arc_df$width_val <- rescale(arc_df$coaccess, to = c(0.8, 2.8))
  } else {
    arc_df$curvature <- 0.45; arc_df$alpha_val <- 0.7; arc_df$width_val <- 1.5
  }

  # --- Prepare Genomic Window ---
  ccan_region <- subset(region_df, ccan == ccan_id)
  if(nrow(ccan_region) == 0) return(NULL)

  region_chr   <- ccan_region$chromosome
  # Namespaced by dir_prefix (global, set in Section 0's CLI parsing): CCAN
  # ids are just local Louvain community labels, not stable across runs (see
  # GRN/CLAUDE.md), so two runs' plots for "CCAN 42" are unrelated regions --
  # writing both to the same un-prefixed 4_figures/<chr>/CCAN_42_track_plot.pdf
  # would silently overwrite one with the other. Matters most for multiple
  # concurrent per-celltype runs (see run_pipeline.sh --celltype-array).
  chr_folder <- file.path(fig_dir, dir_prefix, region_chr)
  if(!dir.exists(chr_folder)) dir.create(chr_folder, recursive = TRUE)

  region_start <- ccan_region$start - 5000
  region_end   <- ccan_region$end + 5000

  # --- PLOT 1: Peaks ---
  peaks_split_all <- do.call(rbind, strsplit(current_peaks, "-"))
  peaks_df <- data.frame(
    chr   = peaks_split_all[, 1],
    start = as.numeric(peaks_split_all[, 2]),
    end   = as.numeric(peaks_split_all[, 3])
  )
  peaks_df$is_hub <- current_peaks %in% top10$Peak

  peak_plot <- ggplot(peaks_df) +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0.2, ymax = 0.8, fill = is_hub)) +
    scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), guide = "none") +
    coord_cartesian(xlim = c(region_start, region_end)) +
    labs(title = paste("ATAC peaks (red = top hubs for CCAN", ccan_id, ")")) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 2)))

  # --- PLOT 2: Arcs ---
  arc_plot <- ggplot() + coord_cartesian(xlim = c(region_start, region_end), ylim = c(0, 1))

  for (i in seq_len(nrow(arc_df))) {
    arc_plot <- arc_plot +
      geom_curve(
        data = arc_df[i, ],
        aes(x = x1, y = 0, xend = x2, yend = 0),
        curvature  = -arc_df$curvature[i],
        alpha      = arc_df$alpha_val[i],
        linewidth  = arc_df$width_val[i],
        color      = "steelblue4"
      )
  }

  arc_plot <- arc_plot +
    labs(title = paste0("CCAN ", ccan_id, " Co-Accessibility (>= 0.95 quant: ", round(threshold, 3), ")")) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 2)))

  # --- PLOT 3: Genes ---
  genes_in_region <- subset(
    gene_gr,
    as.character(seqnames(gene_gr)) == region_chr &
      start(gene_gr) < region_end &
      end(gene_gr)   > region_start
  )

  genes_df <- as.data.frame(genes_in_region)

  if(nrow(genes_df) > 0) {
    genes_df$is_hub_gene <- genes_df$gene_name %in% hub_genes_to_highlight
    genes_df$gene_name <- factor(genes_df$gene_name, levels = unique(genes_df$gene_name[order(genes_df$start)]))
    genes_plus  <- subset(genes_df, strand == "+")
    genes_minus <- subset(genes_df, strand == "-")

    gene_plot <- ggplot() +
      { if (nrow(genes_plus) > 0)
        geom_segment(data = genes_plus, aes(x = start, xend = end, y = gene_name, yend = gene_name, color = is_hub_gene),
                     linewidth = 3, arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "open"))
      } +
      { if (nrow(genes_minus) > 0)
        geom_segment(data = genes_minus, aes(x = start, xend = end, y = gene_name, yend = gene_name, color = is_hub_gene),
                     linewidth = 3, arrow = arrow(length = unit(0.2, "cm"), ends = "first", type = "closed"))
      } +
      scale_color_manual(values = c(`FALSE` = "grey40", `TRUE` = "firebrick"), guide = "none") +
      coord_cartesian(xlim = c(region_start, region_end)) +
      labs(title = "Genes in Region (red = nearest to hubs)", x = paste0("Position: ", region_chr, " (bp)")) +
      theme_minimal(base_size = 9) +
      theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 2)),
            axis.title.y = element_blank())
  } else {
    # Empty placeholder if no genes fall in this window
    gene_plot <- ggplot() + theme_void() + labs(title = "No genes mapped in this region")
  }

  # --- COMBINE AND SAVE ---

  num_genes <- nrow(genes_df)

  dynamic_height <- 5 + (num_genes * 0.25)
  dynamic_height <- max(7, dynamic_height)

  gene_panel_weight <- max(1.8, num_genes * 0.12)

  final_plot <- arc_plot / peak_plot / gene_plot +
    plot_layout(heights = c(2, 0.6, gene_panel_weight))

  file_name <- file.path(chr_folder, paste0("CCAN_", ccan_id, "_track_plot.pdf"))

  ggsave(
    filename = file_name,
    plot = final_plot,
    width = 11,
    height = dynamic_height,
    limitsize = FALSE
  )

  message(paste("Saved PDF plot for CCAN:", ccan_id, "| Genes:", num_genes, "| Height:", round(dynamic_height, 1)))
}
```

### 6. Execute Loop Across All CCANs

The payoff for generalizing Section 5: now plotting *every* CCAN in the dataset is just a `for` loop that calls `plot_ccan_track()` once per CCAN ID. The only thing worth calling out is the `tryCatch()` wrapper — with hundreds/thousands of CCANs, some are going to hit an edge case (weird coordinates, no genes nearby, whatever), and I don't want one bad CCAN to kill a run that's already halfway through everything else. So instead of crashing, a failed CCAN just prints an error message and the loop moves on to the next one.

```r
# Create figure directory if it doesn't exist
if(!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

all_ccan_ids <- unique(ccans$CCAN)
total_ccans <- length(all_ccan_ids)

cat("Starting to generate plots for", total_ccans, "CCANs...\n")

for (i in 1:total_ccans) {
  current_id <- all_ccan_ids[i]

  # Wrap in tryCatch so one bad CCAN doesn't break the whole loop
  tryCatch({
    plot_ccan_track(
      ccan_id   = current_id,
      ccans_df  = ccans,
      conns_df  = conns,
      region_df = region_df,
      gene_gr   = gene_gr,
      fig_dir   = figure_dir
    )
  }, error = function(e) {
    message(paste("Error processing CCAN", current_id, ":", e$message))
  })
}

cat("Finished processing all CCANs.\n")

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
```

I renamed the barcodes to be more easy to identify them using the pattern that pycistopic expects.

```R
cat("=== EXPORTING BARCODES ===\n")

valid_barcodes <- colnames(seurat_obj)
renamed_barcodes <- strsplit(valid_barcodes, split = "_", fixed = TRUE)
ncol(seurat_obj)

for(i in seq_len(ncol(seurat_obj))) {
  print(paste("Processing barcode", i, "of", ncol(seurat_obj))  )
  dataset_name <- seurat_obj$sample_id[i]
  barcode <- renamed_barcodes[[i]][1]
  new_barcode <- paste0(barcode, "_", dataset_name)
  valid_barcodes[i] <- new_barcode
}

colnames(seurat_obj) <- valid_barcodes
writeLines(valid_barcodes, file.path(outdir, "valid_barcodes.txt"))
cat("OK: Exported", length(valid_barcodes), "barcodes\n\n")
```

### Step 3: Examine Metadata

This step is just to verify the data on our object

```R
cat("=== EXAMINING METADATA ===\n")
cat("Columns in metadata:\n")
print(colnames(seurat_obj@meta.data))
cat("\n")

cat("Unique datasets:\n")
print(table(seurat_obj$dataset))
cat("\n")
```

### Step 4: Extract and Verify Fragment Paths

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
cat("=== CREATING SAMPLE MAPPING ===\n")

# Get unique datasets from metadata
unique_datasets <- unique(seurat_obj$dataset)
unique_datasets <- unique_datasets[!is.na(unique_datasets)]
unique_datasets <- sort(unique_datasets)

cat("Unique datasets found:", length(unique_datasets), "\n")
print(unique_datasets)
cat("\n")
```

**Create mapping: dataset name -> fragment path**

IMPORTANT: do NOT assume sort(unique(seurat_obj$dataset)) is in the same order as Fragments(seurat_obj[["ATAC"]]) (frag_list/frag_paths) -- it is not guaranteed to be, and for this object it silently wasn't: every one of the 21 primary ATLASxx_atac samples ended up paired with a DIFFERENT sample's fragments file (confirmed 2026-08-04 via 1_scripts/general/check_fragment_order.r, which reads the ground-truth correspondence directly from the source Seurat object. Instead, determine each fragment file's true dataset directly: each Fragment object's own `cells` slot records which of THIS Seurat object's cells actually draw from it (set when the object was created) -- that's authoritative, not something to infer from list order.

```R
cat("Determining true dataset per fragment file from each Fragment object's own `cells` slot...\n")
fragment_true_dataset <- vapply(frag_list, function(fobj) {
  frag_cells <- intersect(names(fobj@cells), colnames(seurat_obj))
  if (length(frag_cells) == 0) return(NA_character_)
  ds_tab <- sort(table(seurat_obj$dataset[frag_cells]), decreasing = TRUE)
  if (length(ds_tab) > 1 && ds_tab[1] < 0.95 * sum(ds_tab)) {
    stop(sprintf(
      "Fragment file's cells don't cleanly belong to one dataset (top: %s = %d/%d) -- investigate before proceeding.",
      names(ds_tab)[1], ds_tab[1], sum(ds_tab)
    ))
  }
  names(ds_tab)[1]
}, character(1))

if (anyNA(fragment_true_dataset)) {
  stop("Some fragment files had zero cells overlapping seurat_obj -- investigate before proceeding.")
}
if (length(unique(fragment_true_dataset)) != length(frag_paths)) {
  stop("Derived dataset labels for fragment files are not all unique -- investigate before proceeding.")
}

mapping_df <- data.frame(
  sample_id = seq_along(frag_paths),
  dataset = fragment_true_dataset,
  fragment_path = frag_paths,
  stringsAsFactors = FALSE
)

cat("Sample mapping:\n")
print(mapping_df)
cat("\n")

write.table(mapping_df, file.path(outdir, "sample_fragment_mapping.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("OK: Mapping saved to sample_fragment_mapping.tsv\n\n")

```

### Step 6: Stick the Fragment Path onto Every Cell (and Export)

Ok, last one, I promise. We already built `mapping_df`, our little "who owns which fragments file" dictionary. Now we just glue that onto the metadata: for every single cell, look up its `dataset`, find that dataset's row in `mapping_df`, and copy the `fragment_path` over. That's it, that's the whole idea — `match()` does the lookup, so a one-liner does what would otherwise be another one of those `for` loops.

```R
cat("=== ADDING SAMPLE INFO TO METADATA ===\n")

# Map dataset to fragment_path using the mapping_df
seurat_obj$fragment_path <- mapping_df$fragment_path[match(seurat_obj$dataset, mapping_df$dataset)]

cat("Sample assignment summary:\n")
print(table(seurat_obj$dataset))
cat("\n")

cat("Fragment path assignment check:\n")
print(table(!is.na(seurat_obj$fragment_path)))
cat("\n")
```

That last `table(!is.na(...))` is just a cheap sanity check: if any cell came back `FALSE` (i.e. `NA`), it means that cell's dataset had no match in `mapping_df` — which shouldn't happen, but I'd rather see it printed loud and clear than find out three steps later when pycisTopic silently drops cells.

(Optional) Since I only care about Monocytes and Dendritic cells for this GRN, this is also a convenient spot to subset the object down to just those, before writing anything out. Skip this chunk if you want the full `predicted.l1` population instead.

```R
# =========================================================
#  7.6 SUBSET FOR MONOCYTES and DENDRITIC CELLS (OPTIONAL) 
# =========================================================

seurat_obj <- subset(seurat_obj, subset = predicted.l1 %in% c("Monocytes", "Dendritic cells"))
```

And finally, the payoff: dump the whole metadata table (now carrying its brand-new `fragment_path` column) to disk as a tsv. This is deliverable #4, `seurat_metadata_with_fragments.tsv`, the last of the four chivas from the table at the top of this section.

```R
cat("=== EXPORTING METADATA ===\n")

# Full metadata
write.table(seurat_obj@meta.data, file.path(outdir, "seurat_metadata_with_fragments.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
cat("OK: Full metadata exported\n")
```

And with that, all four deliverables are sitting in `2_data/`: `integrated_peaks.bed`, `valid_barcodes.txt`, `sample_fragment_mapping.tsv`, and `seurat_metadata_with_fragments.tsv`. That's everything pycisTopic needs to get going.

## FUTURE WORK

Now I will try to complete pycistopic tutorial. At this point if you have a seurat object with the previous data you could follow it too. btw cloud computing sucks buy me a home lab plspls #grant #jovenpromesa