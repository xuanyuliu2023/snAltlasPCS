# Seurat integration and downstream analysis template
# Replace the placeholder values in the configuration section before running.

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(harmony)
  library(ape)
})

# Configuration
input_object_file <- "INPUT_OBJECT_FILE"
keep_cells_file <- "KEEP_CELLS_FILE"
keep_cells_object_name <- "KEEP_CELL_IDS_OBJECT"
sample_ids <- c("SAMPLE_ID_1", "SAMPLE_ID_2")
output_directory <- "OUTPUT_DIRECTORY"

metadata_file <- ""
gene_type_file <- ""
# Optional: set named paths to pathway matrices to add them as assays.
pathway_files <- character()

normalization_method <- "SCT"
regression_variables <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mito",
  "S.Score",
  "G2M.Score"
)
cluster_resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.6)
marker_resolution <- 0.6
conserved_marker_identity <- ""
conserved_marker_column <- "lineage"
conserved_marker_group <- "group"
discarded_cluster_ids <- character()

max_future_globals_size <- 20 * 1024^3
parallel_workers <- 10

# Validate configuration and create the output directory.
required_files <- c(input_object_file, keep_cells_file)
if (any(grepl("^(INPUT_OBJECT_FILE|KEEP_CELLS_FILE)$", required_files))) {
  stop("Set input_object_file and keep_cells_file before running.")
}
if (length(sample_ids) < 2L ||
    any(grepl("^SAMPLE_ID_", sample_ids)) ||
    anyDuplicated(sample_ids) > 0L) {
  stop("Set sample_ids to the sample identifiers in the Seurat object list.")
}
if (keep_cells_object_name == "KEEP_CELL_IDS_OBJECT") {
  stop("Set keep_cells_object_name to the object name in keep_cells_file.")
}
if (identical(output_directory, "OUTPUT_DIRECTORY")) {
  stop("Set output_directory before running.")
}
if (!normalization_method %in% c("SCT", "LogNormalize")) {
  stop("normalization_method must be either 'SCT' or 'LogNormalize'.")
}
if (length(cluster_resolutions) == 0L ||
    !marker_resolution %in% cluster_resolutions) {
  stop("marker_resolution must be included in cluster_resolutions.")
}
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}
if (!dir.exists(output_directory)) {
  stop("Could not create output_directory.")
}
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Input file(s) not found: ", paste(missing_files, collapse = ", "))
}

plan("multisession", workers = parallel_workers)
options(future.globals.maxSize = max_future_globals_size)

# Load the Seurat object list and the cell identifiers to retain.
input_environment <- new.env(parent = emptyenv())
loaded_object_names <- load(input_object_file, envir = input_environment)
is_seurat_object_list <- vapply(loaded_object_names, function(object_name) {
  candidate <- get(object_name, envir = input_environment)
  is.list(candidate) &&
    length(candidate) > 0L &&
    all(vapply(candidate, inherits, logical(1), what = "Seurat"))
}, logical(1))

seurat_list_names <- loaded_object_names[is_seurat_object_list]
if (length(seurat_list_names) != 1L) {
  stop("The input file must contain exactly one list of Seurat objects.")
}
seurat_list <- get(seurat_list_names, envir = input_environment)

keep_environment <- new.env(parent = emptyenv())
keep_object_names <- load(keep_cells_file, envir = keep_environment)
if (!keep_cells_object_name %in% keep_object_names) {
  stop("keep_cells_object_name was not found in keep_cells_file.")
}
keep_cell_ids <- get(keep_cells_object_name, envir = keep_environment)
if (!is.character(keep_cell_ids)) {
  stop("The configured cell identifier object must be a character vector.")
}

missing_samples <- setdiff(sample_ids, names(seurat_list))
if (length(missing_samples) > 0L) {
  stop("Sample identifiers not found in the Seurat object list: ",
       paste(missing_samples, collapse = ", "))
}
seurat_list <- seurat_list[sample_ids]

# Keep only configured cells in each sample and merge the objects.
seurat_list <- lapply(names(seurat_list), function(sample_name) {
  seurat_object <- seurat_list[[sample_name]]
  cell_ids <- intersect(colnames(seurat_object), keep_cell_ids)
  if (length(cell_ids) == 0L) {
    stop("At least one selected sample has no cells matching keep_cell_ids.")
  }
  seurat_object$sample_id <- sample_name
  subset(seurat_object, cells = cell_ids)
})
names(seurat_list) <- sample_ids

seurat_integrated <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  merge.data = TRUE,
  project = "IntegratedProject"
)

# Normalize the merged object using exactly one configured workflow.
DefaultAssay(seurat_integrated) <- "RNA"
if (normalization_method == "SCT") {
  seurat_integrated <- NormalizeData(
    seurat_integrated,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  seurat_integrated <- CellCycleScoring(
    seurat_integrated,
    s.features = cc.genes$s.genes,
    g2m.features = cc.genes$g2m.genes,
    set.ident = FALSE
  )
  missing_regression_variables <- setdiff(
    regression_variables,
    colnames(seurat_integrated[[]])
  )
  if (length(missing_regression_variables) > 0L) {
    stop("Regression metadata is missing: ",
         paste(missing_regression_variables, collapse = ", "))
  }
  seurat_integrated <- SCTransform(
    seurat_integrated,
    vars.to.regress = regression_variables,
    verbose = FALSE
  )
  integration_assay <- "SCT"
} else {
  seurat_integrated <- NormalizeData(
    seurat_integrated,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  seurat_integrated <- FindVariableFeatures(
    seurat_integrated,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  seurat_integrated <- CellCycleScoring(
    seurat_integrated,
    s.features = cc.genes$s.genes,
    g2m.features = cc.genes$g2m.genes,
    set.ident = FALSE
  )
  missing_regression_variables <- setdiff(
    regression_variables,
    colnames(seurat_integrated[[]])
  )
  if (length(missing_regression_variables) > 0L) {
    stop("Regression metadata is missing: ",
         paste(missing_regression_variables, collapse = ", "))
  }
  seurat_integrated <- ScaleData(
    seurat_integrated,
    vars.to.regress = regression_variables,
    model.use = "linear",
    do.scale = TRUE,
    do.center = TRUE,
    verbose = FALSE
  )
  integration_assay <- "RNA"
}

# Run dimensionality reduction, Harmony integration, and clustering.
seurat_integrated <- RunPCA(
  seurat_integrated,
  assay = integration_assay,
  npcs = 50,
  verbose = FALSE
)
if (!"sample_id" %in% colnames(seurat_integrated[[]])) {
  stop("The metadata must contain a 'sample_id' column for Harmony.")
}
seurat_integrated <- RunHarmony(
  seurat_integrated,
  group.by.vars = "sample_id",
  assay.use = integration_assay
)
seurat_integrated <- RunUMAP(
  seurat_integrated,
  reduction = "harmony",
  dims = 1:30
)
seurat_integrated <- RunTSNE(
  seurat_integrated,
  reduction = "harmony",
  dims = 1:30
)
seurat_integrated <- FindNeighbors(
  seurat_integrated,
  reduction = "harmony",
  dims = 1:30,
  graph.name = c("harmony_nn", "harmony_snn")
)
for (resolution in cluster_resolutions) {
  seurat_integrated <- FindClusters(
    seurat_integrated,
    resolution = resolution,
    graph.name = "harmony_snn"
  )
}

# Add sample-level metadata when a metadata file is configured.
if (nzchar(metadata_file)) {
  if (!file.exists(metadata_file)) {
    stop("metadata_file was not found: ", metadata_file)
  }
  sample_metadata <- read.table(
    metadata_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_metadata_columns <- c("dataID", "category2", "sex", "age")
  missing_metadata_columns <- setdiff(
    required_metadata_columns,
    colnames(sample_metadata)
  )
  if (length(missing_metadata_columns) > 0L) {
    stop("Metadata columns are missing: ",
         paste(missing_metadata_columns, collapse = ", "))
  }
  metadata_index <- match(
    seurat_integrated$sample_id,
    sample_metadata$dataID
  )
  if (anyNA(metadata_index)) {
    stop("Some sample identifiers are missing from metadata_file.")
  }
  seurat_integrated$category <- sample_metadata$category2[metadata_index]
  seurat_integrated$sex <- sample_metadata$sex[metadata_index]
  seurat_integrated$age <- sample_metadata$age[metadata_index]
}

# Find positive markers for the selected clustering resolution.
marker_identity_column <- paste0("harmony_snn_res.", marker_resolution)
if (!marker_identity_column %in% colnames(seurat_integrated[[]])) {
  stop("The selected marker_resolution was not generated.")
}
Idents(seurat_integrated) <- marker_identity_column
markers_all <- FindAllMarkers(
  object = seurat_integrated,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "bimod"
)
write.table(
  markers_all,
  file = file.path(output_directory, "markers_all.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  col.names = TRUE
)

# Retain protein-coding markers when a gene-type annotation file is provided.
if (nzchar(gene_type_file)) {
  if (!file.exists(gene_type_file)) {
    stop("gene_type_file was not found: ", gene_type_file)
  }
  gene_type_info <- read.table(
    gene_type_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  if (ncol(gene_type_info) < 2L) {
    stop("gene_type_file must contain gene identifiers and gene types.")
  }
  protein_coding_genes <- gene_type_info[
    gene_type_info[[2]] == "protein_coding",
    1
  ]
  markers_protein_coding <- subset(
    markers_all,
    gene %in% protein_coding_genes
  )
  write.table(
    markers_protein_coding,
    file = file.path(output_directory, "markers_protein_coding.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    col.names = TRUE
  )
}

# Find conserved markers when both required metadata fields are available.
if (nzchar(conserved_marker_identity)) {
  if (!all(c(conserved_marker_column, conserved_marker_group) %in%
           colnames(seurat_integrated[[]]))) {
    stop("Conserved-marker metadata columns are missing.")
  }
  Idents(seurat_integrated) <- conserved_marker_column
  conserved_markers <- FindConservedMarkers(
    seurat_integrated,
    assay = integration_assay,
    ident.1 = conserved_marker_identity,
    grouping.var = conserved_marker_group,
    verbose = FALSE
  )
  write.table(
    conserved_markers,
    file = file.path(output_directory, "conserved_markers.tsv"),
    sep = "\t",
    row.names = TRUE,
    quote = FALSE,
    col.names = NA
  )
}

# Build and export the cluster tree.
seurat_integrated <- BuildClusterTree(
  object = seurat_integrated,
  assay = "RNA",
  features = VariableFeatures(seurat_integrated, assay = "RNA")
)
cluster_tree <- Tool(seurat_integrated, slot = "BuildClusterTree")
write.tree(
  phy = cluster_tree,
  file = file.path(output_directory, "cluster_tree.newick")
)

# Export expression data for downstream signature scoring.
expression_matrix <- t(as.matrix(
  GetAssayData(seurat_integrated, assay = "RNA", slot = "data")
))
write.table(
  expression_matrix,
  file = file.path(output_directory, "expression_matrix.tsv"),
  sep = "\t",
  row.names = TRUE,
  col.names = NA,
  quote = FALSE
)

# Add pathway score assays when pathway files are configured.
if (length(pathway_files) > 0L) {
  if (is.null(names(pathway_files)) ||
      any(names(pathway_files) == "")) {
    stop("pathway_files must be a named character vector of file paths.")
  }
  for (pathway_name in names(pathway_files)) {
    pathway_file <- pathway_files[[pathway_name]]
    if (!file.exists(pathway_file)) {
      stop("Pathway file was not found: ", pathway_file)
    }
    pathway_matrix <- as.matrix(read.table(
      pathway_file,
      sep = "\t",
      header = TRUE,
      row.names = 1,
      check.names = FALSE
    ))
    pathway_matrix <- t(pathway_matrix)
    if (!setequal(colnames(pathway_matrix), colnames(seurat_integrated))) {
      stop("Pathway matrix cell identifiers do not match the Seurat object.")
    }
    pathway_matrix <- pathway_matrix[, colnames(seurat_integrated), drop = FALSE]
    assay_name <- paste0("PATHWAY_", toupper(pathway_name))
    seurat_integrated[[assay_name]] <- CreateAssayObject(
      counts = pathway_matrix
    )
  }
}

# Save the integrated object and optional cell lists.
save(
  seurat_integrated,
  file = file.path(output_directory, "seurat_integrated.RData")
)

if (length(discarded_cluster_ids) > 0L) {
  discard_column <- "harmony_snn_res.1"
  if (!discard_column %in% colnames(seurat_integrated[[]])) {
    stop("The configured discard clusters require resolution 1.0.")
  }
  cells_to_discard <- rownames(seurat_integrated[[]])[
    seurat_integrated[[discard_column, drop = TRUE]] %in%
      discarded_cluster_ids
  ]
  save(
    cells_to_discard,
    file = file.path(output_directory, "cells_to_discard.RData")
  )
}

cells_to_keep <- rownames(seurat_integrated[[]])[
  seurat_integrated$sample_id %in% sample_ids
]
save(
  cells_to_keep,
  file = file.path(output_directory, "cells_to_keep.RData")
)
