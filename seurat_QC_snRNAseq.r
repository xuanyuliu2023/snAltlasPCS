# Single-cell RNA-seq quality control template
# Replace the placeholder values in the configuration section before running.

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(Matrix)
  library(ggplot2)
  library(ggpubr)
})

# Configuration
species <- "human"
sample_ids <- c("SAMPLE_ID_1", "SAMPLE_ID_2")
cellranger_root <- "CELLRANGER_ROOT"
cellbender_root <- "CELLBENDER_ROOT"
output_directory <- "OUTPUT_DIRECTORY"

cellranger_relative_path <- file.path(
  "{sample_id}", "outs", "filtered_feature_bc_matrix.h5"
)
cellbender_relative_path <- file.path(
  "{sample_id}", "{sample_id}_raw_feature_bc_matrix_cellbender_filtered_seurat.h5"
)

# Set a gene-list file to use a custom ribosomal gene list; otherwise use a prefix.
ribosomal_gene_file <- ""
ribosomal_gene_prefix <- if (species == "human") "^(RPL|RPS)" else "^(Rpl|Rps)"

# Cell-cycle scoring is optional. The RData file must contain one object with
# character vectors named s.genes and g2m.genes.
cell_cycle_gene_file <- ""

# Set a directory to remove predicted doublets after exporting matrix files.
# Each prediction file must contain one comma-separated TRUE/FALSE value per cell.
doublet_prediction_directory <- ""

parallel_workers <- 10
max_future_globals_size <- 20 * 1024^3

# Enter numeric cutoffs before running. Values apply to every selected sample.
qc_thresholds <- data.frame(
  metric = c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo"),
  lower = c(NA_real_, NA_real_, NA_real_, NA_real_),
  upper = c(NA_real_, NA_real_, NA_real_, NA_real_),
  stringsAsFactors = FALSE
)

# Validate configuration and required inputs.
if (!species %in% c("human", "mouse")) {
  stop("species must be either 'human' or 'mouse'.")
}
if (length(sample_ids) == 0L ||
    any(grepl("^SAMPLE_ID_", sample_ids)) ||
    anyDuplicated(sample_ids) > 0L) {
  stop("Set sample_ids to unique sample identifiers.")
}
if (any(grepl("^(CELLRANGER_ROOT|CELLBENDER_ROOT|OUTPUT_DIRECTORY)$",
              c(cellranger_root, cellbender_root, output_directory)))) {
  stop("Set cellranger_root, cellbender_root, and output_directory.")
}
required_qc_metrics <- c(
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mito",
  "percent.ribo"
)
if (anyDuplicated(qc_thresholds$metric) > 0L ||
    !setequal(qc_thresholds$metric, required_qc_metrics)) {
  stop("qc_thresholds must define each required metric exactly once.")
}
if (anyNA(qc_thresholds$lower) ||
    anyNA(qc_thresholds$upper) ||
    any(qc_thresholds$lower >= qc_thresholds$upper)) {
  stop("Set valid lower and upper values for every QC metric.")
}

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}
if (!dir.exists(output_directory)) {
  stop("Could not create output_directory.")
}

plan("multisession", workers = parallel_workers)
options(future.globals.maxSize = max_future_globals_size)

# Read a single gene-expression matrix from a 10x HDF5 file.
read_expression_matrix <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("Input file was not found: ", file_path)
  }
  matrix_data <- Read10X_h5(file_path)
  if (is.list(matrix_data)) {
    if ("Gene Expression" %in% names(matrix_data)) {
      matrix_data <- matrix_data[["Gene Expression"]]
    } else if (length(matrix_data) == 1L) {
      matrix_data <- matrix_data[[1]]
    } else {
      stop("Could not identify a unique gene-expression matrix in: ",
           file_path)
    }
  }
  if (is.null(rownames(matrix_data)) || is.null(colnames(matrix_data))) {
    stop("The input matrix must have gene and cell identifiers.")
  }
  matrix_data
}

# Add a sample identifier to every cell barcode.
make_cell_names <- function(cell_ids, sample_id) {
  paste(cell_ids, sample_id, sep = "___")
}

# Import counts, retain shared barcodes, and create one object per sample.
seurat_objects <- setNames(vector("list", length(sample_ids)), sample_ids)
for (sample_id in sample_ids) {
  cellranger_file <- file.path(
    cellranger_root,
    gsub("\\{sample_id\\}", sample_id, cellranger_relative_path)
  )
  cellbender_file <- file.path(
    cellbender_root,
    gsub("\\{sample_id\\}", sample_id, cellbender_relative_path)
  )

  cellranger_counts <- read_expression_matrix(cellranger_file)
  cellbender_counts <- read_expression_matrix(cellbender_file)
  shared_cells <- intersect(
    colnames(cellranger_counts),
    colnames(cellbender_counts)
  )
  if (length(shared_cells) == 0L) {
    stop("No shared cell barcodes were found for sample: ", sample_id)
  }

  cellranger_counts <- cellranger_counts[, shared_cells, drop = FALSE]
  cellbender_counts <- cellbender_counts[, shared_cells, drop = FALSE]
  renamed_cells <- make_cell_names(shared_cells, sample_id)
  colnames(cellranger_counts) <- renamed_cells
  colnames(cellbender_counts) <- renamed_cells

  seurat_object <- CreateSeuratObject(
    counts = cellranger_counts,
    project = sample_id,
    assay = "RNA_cellranger",
    min.cells = 3,
    min.features = 0
  )
  seurat_object[["RNA_cellbender"]] <- CreateAssayObject(
    counts = cellbender_counts
  )
  seurat_object$sample_id <- sample_id

  # Calculate QC counts from the assay used for downstream analysis.
  cellbender_counts <- GetAssayData(
    seurat_object,
    assay = "RNA_cellbender",
    slot = "counts"
  )
  seurat_object$nCount_RNA <- Matrix::colSums(cellbender_counts)
  seurat_object$nFeature_RNA <- Matrix::colSums(cellbender_counts > 0)
  DefaultAssay(seurat_object) <- "RNA_cellbender"

  if (species == "human") {
    mitochondrial_pattern <- "^MT-"
    rbc_pattern <- "^HB[ABDEGMQZ]"
  } else {
    mitochondrial_pattern <- "^mt-"
    rbc_pattern <- "^Hb[abdeqxyz]"
  }
  seurat_object[["percent.mito"]] <- PercentageFeatureSet(
    seurat_object,
    pattern = mitochondrial_pattern,
    assay = "RNA_cellbender"
  )
  seurat_object[["percent.rbc"]] <- PercentageFeatureSet(
    seurat_object,
    pattern = rbc_pattern,
    assay = "RNA_cellbender"
  )

  if (nzchar(ribosomal_gene_file)) {
    if (!file.exists(ribosomal_gene_file)) {
      stop("ribosomal_gene_file was not found: ", ribosomal_gene_file)
    }
    ribosomal_genes <- read.table(
      ribosomal_gene_file,
      header = FALSE,
      stringsAsFactors = FALSE
    )[[1]]
    ribosomal_features <- intersect(
      ribosomal_genes,
      rownames(seurat_object[["RNA_cellbender"]])
    )
    if (length(ribosomal_features) == 0L) {
      stop("No ribosomal genes matched sample: ", sample_id)
    }
    seurat_object[["percent.ribo"]] <- PercentageFeatureSet(
      seurat_object,
      features = ribosomal_features,
      assay = "RNA_cellbender"
    )
  } else {
    seurat_object[["percent.ribo"]] <- PercentageFeatureSet(
      seurat_object,
      pattern = ribosomal_gene_prefix,
      assay = "RNA_cellbender"
    )
  }

  seurat_object <- NormalizeData(
    seurat_object,
    assay = "RNA_cellbender",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  seurat_objects[[sample_id]] <- seurat_object
}

# Optionally score cell-cycle phases from a configured gene-set file.
if (nzchar(cell_cycle_gene_file)) {
  if (!file.exists(cell_cycle_gene_file)) {
    stop("cell_cycle_gene_file was not found: ", cell_cycle_gene_file)
  }
  cell_cycle_environment <- new.env(parent = emptyenv())
  cell_cycle_object_names <- load(
    cell_cycle_gene_file,
    envir = cell_cycle_environment
  )
  cell_cycle_candidates <- cell_cycle_object_names[vapply(
    cell_cycle_object_names,
    function(object_name) {
      candidate <- get(object_name, envir = cell_cycle_environment)
      is.list(candidate) &&
        all(c("s.genes", "g2m.genes") %in% names(candidate))
    },
    logical(1)
  )]
  if (length(cell_cycle_candidates) != 1L) {
    stop("The cell-cycle file must contain exactly one compatible gene set.")
  }
  cell_cycle_genes <- get(
    cell_cycle_candidates,
    envir = cell_cycle_environment
  )
  for (sample_id in sample_ids) {
    seurat_objects[[sample_id]] <- CellCycleScoring(
      seurat_objects[[sample_id]],
      s.features = cell_cycle_genes$s.genes,
      g2m.features = cell_cycle_genes$g2m.genes,
      assay = "RNA_cellbender",
      set.ident = FALSE
    )
  }
}

# Save unfiltered objects and export per-sample QC plots and summaries.
save(
  seurat_objects,
  file = file.path(output_directory, "seurat_objects_unfiltered.RData")
)

qc_summary <- data.frame(
  sample_id = sample_ids,
  cells_before_qc = vapply(
    seurat_objects,
    ncol,
    integer(1)
  ),
  cells_after_qc = integer(length(sample_ids)),
  stringsAsFactors = FALSE
)

make_qc_scatterplot <- function(metadata, metric_name) {
  count_threshold <- qc_thresholds[
    qc_thresholds$metric == "nCount_RNA",
    ,
    drop = FALSE
  ]
  plot <- ggplot(
    metadata,
    aes(x = nCount_RNA, y = .data[[metric_name]])
  ) +
    geom_point(size = 0.7, alpha = 0.7, color = "#2E9FDF") +
    geom_vline(
      xintercept = count_threshold$lower,
      linetype = "dashed",
      color = "red"
    ) +
    geom_vline(
      xintercept = count_threshold$upper,
      linetype = "dashed",
      color = "red"
    ) +
    labs(
      title = metric_name,
      x = "Total RNA counts",
      y = metric_name
    ) +
    theme_bw()

  metric_threshold <- qc_thresholds[
    qc_thresholds$metric == metric_name,
    ,
    drop = FALSE
  ]
  if (nrow(metric_threshold) == 1L) {
    plot <- plot +
      geom_hline(
        yintercept = metric_threshold$lower,
        linetype = "dashed",
        color = "red"
      ) +
      geom_hline(
        yintercept = metric_threshold$upper,
        linetype = "dashed",
        color = "red"
      )
  }
  plot
}

for (sample_id in sample_ids) {
  sample_output_directory <- file.path(output_directory, sample_id)
  if (!dir.exists(sample_output_directory)) {
    dir.create(sample_output_directory, recursive = TRUE)
  }

  metadata <- seurat_objects[[sample_id]][[]]
  required_metrics <- qc_thresholds$metric
  if (!all(required_metrics %in% colnames(metadata))) {
    stop("One or more QC metrics are missing for sample: ", sample_id)
  }

  density_plots <- lapply(required_metrics, function(metric_name) {
    ggplot(metadata, aes(x = .data[[metric_name]])) +
      geom_density(fill = "#2E9FDF", alpha = 0.5) +
      labs(title = metric_name, x = metric_name, y = "Density") +
      theme_bw()
  })
  density_figure <- ggarrange(
    plotlist = density_plots,
    ncol = 2,
    nrow = 2
  )
  ggsave(
    filename = file.path(sample_output_directory, "qc_metric_densities.pdf"),
    plot = density_figure,
    width = 10,
    height = 7
  )

  scatter_plots <- lapply(
    c("percent.mito", "nFeature_RNA", "percent.ribo", "percent.rbc"),
    function(metric_name) {
      make_qc_scatterplot(metadata, metric_name)
    }
  )
  scatter_figure <- ggarrange(
    plotlist = scatter_plots,
    ncol = 2,
    nrow = 2
  )
  ggsave(
    filename = file.path(sample_output_directory, "qc_metric_scatterplots.pdf"),
    plot = annotate_figure(
      scatter_figure,
      top = text_grob(sample_id)
    ),
    width = 12,
    height = 9
  )

  write.table(
    metadata,
    file = file.path(sample_output_directory, "cell_metadata.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

# Apply the configured thresholds to each sample.
seurat_objects_qc_filtered <- setNames(
  vector("list", length(sample_ids)),
  sample_ids
)
for (sample_id in sample_ids) {
  metadata <- seurat_objects[[sample_id]][[]]
  keep_cells <- rep(TRUE, nrow(metadata))
  for (row_index in seq_len(nrow(qc_thresholds))) {
    metric_name <- qc_thresholds$metric[row_index]
    keep_cells <- keep_cells &
      metadata[[metric_name]] >= qc_thresholds$lower[row_index] &
      metadata[[metric_name]] <= qc_thresholds$upper[row_index]
  }
  cell_names <- rownames(metadata)[keep_cells]
  if (length(cell_names) == 0L) {
    stop("QC thresholds retained no cells for sample: ", sample_id)
  }
  seurat_objects_qc_filtered[[sample_id]] <- subset(
    seurat_objects[[sample_id]],
    cells = cell_names
  )
  qc_summary$cells_after_qc[
    qc_summary$sample_id == sample_id
  ] <- length(cell_names)
}

# Export filtered matrices and identifiers for an external Scrublet run.
matrix_output_directory <- file.path(output_directory, "scrublet_input")
if (!dir.exists(matrix_output_directory)) {
  dir.create(matrix_output_directory, recursive = TRUE)
}
for (sample_id in sample_ids) {
  count_matrix <- GetAssayData(
    seurat_objects_qc_filtered[[sample_id]],
    assay = "RNA_cellbender",
    slot = "counts"
  )
  writeMM(
    count_matrix,
    file = file.path(matrix_output_directory, paste0(sample_id, ".matrix.mtx"))
  )
  write.table(
    rownames(count_matrix),
    file = file.path(matrix_output_directory, paste0(sample_id, ".genes.tsv")),
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  write.table(
    colnames(count_matrix),
    file = file.path(matrix_output_directory, paste0(sample_id, ".barcodes.tsv")),
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

# Import Scrublet predictions when a prediction directory is configured.
if (nzchar(doublet_prediction_directory)) {
  if (!dir.exists(doublet_prediction_directory)) {
    stop("doublet_prediction_directory does not exist.")
  }
  for (sample_id in sample_ids) {
    prediction_file <- file.path(
      doublet_prediction_directory,
      paste0(sample_id, ".predicted_doublets.csv")
    )
    if (!file.exists(prediction_file)) {
      stop("Scrublet prediction file was not found: ", prediction_file)
    }
    doublet_predictions <- scan(
      prediction_file,
      what = character(),
      sep = ",",
      quiet = TRUE
    )
    expected_cells <- colnames(seurat_objects_qc_filtered[[sample_id]])
    if (length(doublet_predictions) != length(expected_cells) ||
        any(!doublet_predictions %in% c("TRUE", "FALSE"))) {
      stop("Scrublet predictions do not match the cells for sample: ",
           sample_id)
    }
    seurat_objects_qc_filtered[[sample_id]]$doublet <-
      doublet_predictions == "TRUE"
    seurat_objects_qc_filtered[[sample_id]] <- subset(
      seurat_objects_qc_filtered[[sample_id]],
      subset = doublet == FALSE
    )
  }
}

qc_summary$cells_after_doublet_filter <- vapply(
  seurat_objects_qc_filtered,
  ncol,
  integer(1)
)
write.table(
  qc_summary,
  file = file.path(output_directory, "qc_cell_counts.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

save(
  seurat_objects_qc_filtered,
  file = file.path(output_directory, "seurat_objects_qc_filtered.RData")
)
