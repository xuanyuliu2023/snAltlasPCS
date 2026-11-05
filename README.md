# snAltlasPCS
# An integrated single-nucleus and spatial transcriptomics atlas defines lineage plasticity and subtype-specific microenvironment of primary cardiac sarcoma

This repository contains a collection of reusable analysis scripts for a multi-omics single-cell and spatial transcriptomics workflow focused on primary cardiac sarcoma. The scripts cover quality control, integration, copy-number inference, velocity analysis, lineage inference, and spatial profiling.

## Overview

The analysis pipeline is organized around the following major tasks:

- Single-nucleus RNA-seq preprocessing and integration
- Quality control and filtering for single-cell transcriptomics
- Copy-number alteration analysis using SCEVAN
- Spatial transcriptomics processing and cell segmentation
- RNA velocity analysis using scVelo
- Lineage and fate modeling with CellRank
- Drug-target activity scoring with Drug2Cell

## Repository contents

### 1. seurat_QC_snRNAseq.r
Purpose:
- Performs per-sample single-nucleus QC and filtering.
- Calculates standard metrics such as library size, detected genes, mitochondrial percentage, and ribosome percentage.
- Exports filtered matrices for downstream doublet prediction and integration.

Typical workflow:
- Read 10x HDF5 expression matrices
- Harmonize sample barcodes and metadata
- Compute QC metrics
- Apply threshold-based filtering
- Export processed matrices and optional doublet predictions

### 2. seurat_integration.r
Purpose:
- Integrates multiple single-nucleus RNA-seq datasets into a unified Seurat object.
- Keeps only preselected cells and harmonizes batch effects across samples.
- Performs normalization, dimensionality reduction, clustering, and label transfer.

Typical workflow:
- Load a list of Seurat objects
- Filter to retained cell barcodes
- Merge samples using sample identifiers
- Normalize and regress technical covariates
- Run PCA and Harmony integration
- Cluster cells and export marker-based outputs

### 3. SCEVAN.r
Purpose:
- Runs copy-number alteration (CNA) analysis for single-nucleus data.
- Supports both single-sample and multi-sample comparison workflows.
- Generates per-sample result objects and summary plots of clonal states.

Typical workflow:
- Build protein-coding count matrices from the integrated Seurat object
- Define matched normal reference cells
- Run SCEVAN single-sample CNA analysis
- Compare post-processed clones across samples
- Export merged results and visualization plots

### 4. Visium_HD_seurat.r
Purpose:
- Processes and visualizes Visium HD spatial transcriptomics data in Seurat.
- Loads spatial image and segmentation information, performs clustering, and overlays cell annotations.
- Enables zoomed tissue region exploration and feature plotting.

Typical workflow:
- Load Visium HD data with spatial images
- Inspect spatial coordinates, tissue geometry, and segmentation structure
- Generate count summaries and QC plots
- Run PCA, clustering, and UMAP on the spatial assay
- Transfer single-cell labels into spatial object and plot cell-type patterns

### 5. scVelo.py
Purpose:
- Combines Loom-formatted spliced/unspliced counts into a unified RNA velocity dataset.
- Normalizes expression, estimates dynamics, and computes latent time or pseudotime trends.
- Produces velocity stream, grid, and fate-related plots.

Typical workflow:
- Merge multiple Loom files
- Standardize and reorder cell identifiers
- Estimate RNA velocity and transition graphs
- Generate velocity embedding and state plots
- Export an AnnData object and ranked gene summaries

### 6. cellrank.py
Purpose:
- Integrates scVelo outputs with CellRank to infer lineage states and fate probabilities.
- identifies macrostates, terminal and initial states, and lineage-specific driver genes.
- Supports optional gene-trend analysis and output export.

Typical workflow:
- Load an annotated AnnData object
- Preprocess and estimate RNA velocity
- Build connectivity and velocity kernels
- Fit GPCCA lineage models
- Compute fate probabilities and driver genes
- Save lineage summaries and selected plots

### 7. Drug2cell.py
Purpose:
- Evaluates drug-target activity using the Drug2Cell framework on single-cell data.
- Scores drug response patterns, ranks features across groups, and exports summary files.
- Provides optional UMAP visualization and category-based filtering.

Typical workflow:
- Load an AnnData object
- Score drug-target activity for cells or groups
- Rank features by group significance
- Export top genes and drug-category summaries
- Visualize the final score structure on UMAP and save results

## Recommended execution order

For a typical project workflow, the order is:

1. seurat_QC_snRNAseq.r
2. seurat_integration.r
3. SCEVAN.r
4. Visium_HD_seurat.r
5. scVelo.py
6. cellrank.py
7. Drug2cell.py

## Notes

- Several scripts are written as reusable templates and require project-specific configuration values before execution.
- File paths, sample identifiers, metadata columns, and output directories must be customized for each dataset.
- The repository is intended for computational analysis and is designed to be adapted to new datasets with minimal code changes.

## License

This project is distributed for research use. Please check the licensing and data-use constraints of any upstream reference datasets, software tools, or institutional resources used in the analysis.
