# Configuration and package setup
suppressPackageStartupMessages({
  library(Seurat)
  library(SCEVAN)
  library(ggalluvial)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(pheatmap)
  library(plyr)
})

sample_id_1 <- "sample_id_1"
sample_id_2 <- "sample_id_2"
sample_id_3 <- "sample_id_3"
sample_id_4 <- "sample_id_4"
sample_id_5 <- "sample_id_5"
sample_id_6 <- "sample_id_6"
sample_id_7 <- "sample_id_7"
sample_id_8 <- "sample_id_8"

config <- list(
  project_root = "<project_root>",
  integrated_object_path = file.path("<project_root>", "<integrated_seurat_object>.RData"),
  gene_type_file = file.path("<project_root>", "<gene_type_table>.tsv"),
  output_dir = file.path("<project_root>", "<output_directory>"),
  sample_ids = c(sample_id_1, sample_id_2, sample_id_3, sample_id_4, sample_id_5, sample_id_6, sample_id_7, sample_id_8),
  normal_cluster_names = c("cluster_10", "cluster_19", "cluster_4", "cluster_14", "cluster_20"),
  analysis_name = "<analysis_name>",
  organism = "human",
  par_cores = 8L,
  clone_plot = TRUE,
  tree_plot = TRUE,
  beta_vega = 0.5,
  single_sample_id = sample_id_8,
  plot_theme = "minimal"
)

validate_config <- function(config) {
  if (!dir.exists(config$output_dir)) {
    dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (!file.exists(config$integrated_object_path)) {
    stop("The integrated Seurat object is missing.")
  }

  if (!file.exists(config$gene_type_file)) {
    stop("The gene annotation table is missing.")
  }

  if (length(config$sample_ids) == 0) {
    stop("No sample IDs were provided.")
  }
}

resolve_cell_id_column <- function(meta_data) {
  preferred_names <- c("cellID", "cell_id", "cell_barcode", "cell_name")
  match_index <- which(preferred_names %in% colnames(meta_data))
  if (length(match_index) > 0) {
    return(preferred_names[match_index[1]])
  }
  candidate_names <- colnames(meta_data)[grepl("cell|barcode|identifier", colnames(meta_data), ignore.case = TRUE)]
  if (length(candidate_names) > 0) {
    return(candidate_names[1])
  }
  stop("No valid cell identifier column was found in the Seurat metadata.")
}

resolve_cluster_column <- function(meta_data) {
  preferred_names <- c("subcluster", "cluster", "seurat_clusters", "celltype", "cluster_id")
  match_index <- which(preferred_names %in% colnames(meta_data))
  if (length(match_index) > 0) {
    return(preferred_names[match_index[1]])
  }
  candidate_names <- colnames(meta_data)[grepl("cluster|subcluster|type", colnames(meta_data), ignore.case = TRUE)]
  if (length(candidate_names) > 0) {
    return(candidate_names[1])
  }
  stop("No valid cluster annotation column was found in the Seurat metadata.")
}

prepare_count_matrix <- function(seurat_object, gene_type_file) {
  all_genes <- rownames(seurat_object)
  gene_type_info <- read.table(
    file = gene_type_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  if (ncol(gene_type_info) < 2) {
    stop("The gene annotation table must contain at least two columns.")
  }

  protein_coding_genes <- subset(gene_type_info, gene_type_info[, 2] == "protein_coding")[[1]]
  stable_genes <- intersect(all_genes, protein_coding_genes)

  if (length(stable_genes) == 0) {
    stop("No protein-coding genes were retained after intersecting with the annotation table.")
  }

  counts_matrix <- seurat_object@assays$RNA@counts
  count_matrix <- t(as.matrix(counts_matrix[stable_genes, , drop = FALSE]))
  return(count_matrix)
}

save_sample_count_matrix <- function(seurat_object, config) {
  count_matrix <- prepare_count_matrix(seurat_object, config$gene_type_file)
  cell_id_column <- resolve_cell_id_column(seurat_object@meta.data)

  for (sample_id in config$sample_ids) {
    cell_ids <- subset(seurat_object@meta.data, orig.ident == sample_id)[[cell_id_column]]
    if (length(cell_ids) == 0) {
      warning(sprintf("No cells were found for sample %s.", sample_id))
      next
    }

    sample_matrix <- count_matrix[, cell_ids, drop = FALSE]
    output_file <- file.path(config$output_dir, sprintf("countMatrix_cells_%s.RData", sample_id))
    save(sample_matrix, file = output_file)
  }
}

save_normal_cells <- function(seurat_object, config) {
  cell_id_column <- resolve_cell_id_column(seurat_object@meta.data)
  cluster_column <- resolve_cluster_column(seurat_object@meta.data)

  for (sample_id in config$sample_ids) {
    normal_cells <- subset(
      seurat_object@meta.data,
      orig.ident == sample_id & get(cluster_column) %in% config$normal_cluster_names
    )[[cell_id_column]]

    if (length(normal_cells) == 0) {
      warning(sprintf("No normal reference cells were found for sample %s.", sample_id))
      next
    }

    output_file <- file.path(config$output_dir, sprintf("Normal_cells_%s.RData", sample_id))
    save(normal_cells, file = output_file)
  }
}

save_all_cells <- function(seurat_object, config) {
  cell_id_column <- resolve_cell_id_column(seurat_object@meta.data)

  for (sample_id in config$sample_ids) {
    all_cells <- subset(seurat_object@meta.data, orig.ident == sample_id)[[cell_id_column]]
    if (length(all_cells) == 0) {
      warning(sprintf("No cells were found for sample %s.", sample_id))
      next
    }

    output_file <- file.path(config$output_dir, sprintf("ALL_cells_%s.RData", sample_id))
    save(all_cells, file = output_file)

    result_vector <- sapply(strsplit(all_cells, "___"), function(x) x[1])
    output_table <- file.path(config$output_dir, sprintf("ALL_cells_%s.tsv", sample_id))
    write.table(result_vector, file = output_table, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

run_single_sample_cna <- function(seurat_object, config, sample_id) {
  count_path <- file.path(config$output_dir, sprintf("countMatrix_cells_%s.RData", sample_id))
  normal_path <- file.path(config$output_dir, sprintf("Normal_cells_%s.RData", sample_id))

  if (!file.exists(count_path) || !file.exists(normal_path)) {
    stop(sprintf("Missing input files for single-sample CNA analysis of %s.", sample_id))
  }

  count_mtx <- get(load(count_path))
  norm_cells <- get(load(normal_path))

  results <- SCEVAN::pipelineCNA(
    count_mtx,
    norm_cell = norm_cells,
    sample = sample_id,
    par_cores = config$par_cores,
    SUBCLONES = TRUE,
    plotTree = config$tree_plot,
    beta_vega = config$beta_vega,
    ClonalCN = TRUE
  )

  output_file <- file.path(config$output_dir, sprintf("results_%s.RData", sample_id))
  save(results, file = output_file)
  return(results)
}

run_multi_sample_cna <- function(config, sample_list) {
  count_list <- list()
  normal_cells <- list()

  for (sample_id in sample_list) {
    count_path <- file.path(config$output_dir, sprintf("countMatrix_cells_%s.RData", sample_id))
    normal_path <- file.path(config$output_dir, sprintf("Normal_cells_%s.RData", sample_id))

    if (!file.exists(count_path) || !file.exists(normal_path)) {
      stop(sprintf("Missing data for multi-sample CNA comparison: %s.", sample_id))
    }

    count_list[[sample_id]] <- get(load(count_path))
    normal_cells[[sample_id]] <- get(load(normal_path))
  }

  results <- SCEVAN::multiSampleComparisonClonalCN(
    count_list,
    analysisName = config$analysis_name,
    listNormCells = normal_cells,
    organism = config$organism,
    par_cores = config$par_cores
  )

  output_file <- file.path(config$output_dir, sprintf("results_%s_multisample.RData", config$analysis_name))
  save(results, file = output_file)
  return(results)
}

build_alluvial_plot <- function(results, config) {
  if (is.null(results$lineage) || is.null(results$class)) {
    stop("The result table must contain lineage and class columns for alluvial plotting.")
  }

  df <- results[, c("lineage", "class")]
  df$Freq <- 1
  df_wide <- df %>%
    group_by(lineage, class) %>%
    summarize(Freq = n(), .groups = "drop")

  plot <- ggplot(
    data = df_wide,
    aes(axis1 = lineage, axis2 = class, y = Freq)
  ) +
    scale_x_discrete(limits = c("lineage", "class"), expand = c(0.2, 0.05)) +
    xlab("subcluster") +
    geom_alluvium(aes(fill = lineage)) +
    geom_stratum(aes(fill = lineage)) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal() +
    theme(legend.position = "none")

  output_file <- file.path(config$output_dir, "alluvial_plot.pdf")
  ggsave(filename = output_file, plot = plot, device = "pdf")
  return(plot)
}

build_proportion_plot <- function(results, config) {
  if (!all(c("subclone", "subcluster") %in% colnames(results))) {
    stop("The result table must contain subclone and subcluster columns for the proportion plot.")
  }

  subset_results <- subset(results, class == "tumor")
  obs_counts <- as.data.frame(as.matrix(table(subset_results[, c("subclone", "subcluster") ])))
  obs_counts <- dcast(obs_counts, subclone ~ subcluster)
  rownames(obs_counts) <- obs_counts$subclone
  obs_counts$subclone <- paste0("sc", obs_counts$subclone)
  rownames(obs_counts) <- obs_counts$subclone
  obs_counts <- obs_counts[, -1, drop = FALSE]

  prop <- apply(obs_counts, 1, function(x) x / sum(x))
  prop_melt_df <- melt(prop)
  prop_melt_df$class <- rep(c("sc1", "sc2", "sc3", "sc4", "sc5"), 5)
  colnames(prop_melt_df) <- c("subcluster", "Prop", "subclone")

  plot <- ggplot(prop_melt_df, aes(x = subclone, y = Prop * 100, fill = subcluster)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    labs(x = "", y = "Relative proportion (%)")

  output_file <- file.path(config$output_dir, "seurat_cluster_prop_eachSample.pdf")
  ggsave(filename = output_file, plot = plot, device = "pdf")
  return(plot)
}

combine_sample_results <- function(config) {
  result_files <- list.files(
    path = config$output_dir,
    pattern = "^results_.*\\.RData$",
    full.names = TRUE
  )

  if (length(result_files) == 0) {
    stop("No per-sample results were found to combine.")
  }

  combined_results <- NULL
  for (result_file in result_files) {
    sample_result <- get(load(result_file))
    combined_results <- rbind(combined_results, sample_result)
  }

  output_file <- file.path(config$output_dir, "SCEVAN_results_allSamps.RData")
  save(combined_results, file = output_file)
  return(combined_results)
}

main <- function() {
  validate_config(config)

  seurat_object <- get(load(config$integrated_object_path))
  if (class(seurat_object)[1] != "Seurat") {
    stop("The loaded object is not a Seurat object.")
  }

  save_sample_count_matrix(seurat_object, config)
  save_normal_cells(seurat_object, config)
  save_all_cells(seurat_object, config)

  for (sample_id in config$sample_ids) {
    run_single_sample_cna(seurat_object, config, sample_id)
  }

  multi_sample_results <- run_multi_sample_cna(config, config$sample_ids)
  if (!is.null(multi_sample_results)) {
    if (!is.null(multi_sample_results$lineage) && !is.null(multi_sample_results$class)) {
      build_alluvial_plot(multi_sample_results, config)
    }
    if (all(c("subclone", "subcluster", "class") %in% colnames(multi_sample_results))) {
      build_proportion_plot(multi_sample_results, config)
    }
  }

  combined_results <- combine_sample_results(config)
  return(combined_results)
}

results_all_samples <- main()

# End of analysis workflow
