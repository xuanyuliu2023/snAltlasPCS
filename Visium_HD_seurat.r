# Reusable Seurat workflow for Visium HD cell segmentation data.
# Update the configuration values before running the analysis.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(future)
})

# Input and output settings
config <- list(
  sample_id = "sample_id",
  spatial_data_dir = file.path("<project_root>", "<spatial_data_directory>"),
  reference_object_file = file.path("<project_root>", "<reference_object_file>.RData"),
  output_dir = file.path("<project_root>", "<output_directory>"),
  image_name = "tissue_hires_image.png",
  bin_size = "polygons",
  assay = "Spatial.Polygons",
  reference_label_columns = c("label_level_1", "label_level_2"),
  reduction_dims = 1:30,
  clustering_dims = 1:30,
  clustering_resolution = 1,
  workers = 4L,
  max_memory_gb = 50,
  zoom_bounds = NULL,
  features_to_plot = character(),
  segmentation_color = "turquoise4"
)

validate_config <- function(config) {
  if (grepl("^<.*>$", config$spatial_data_dir) ||
      grepl("^<.*>$", config$reference_object_file) ||
      grepl("^<.*>$", config$output_dir)) {
    stop("Replace the path placeholders in config before running the workflow.")
  }

  if (!dir.exists(config$spatial_data_dir)) {
    stop("The configured spatial data directory does not exist.")
  }

  if (!file.exists(config$reference_object_file)) {
    stop("The configured reference object file does not exist.")
  }

  if (length(config$sample_id) != 1L || !nzchar(config$sample_id)) {
    stop("Configure exactly one non-empty sample identifier.")
  }

  if (config$workers < 1L || config$max_memory_gb <= 0) {
    stop("Worker count and memory limit must be positive values.")
  }

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
}

load_single_object <- function(file_path) {
  loaded_names <- load(file_path, envir = environment())
  if (length(loaded_names) != 1L) {
    stop("The input RData file must contain exactly one R object.")
  }
  get(loaded_names[[1L]], envir = environment())
}

save_plot <- function(plot, file_name, config, width = 10, height = 8) {
  output_path <- file.path(config$output_dir, file_name)
  ggsave(filename = output_path, plot = plot, width = width, height = height)
}

load_spatial_object <- function(config) {
  spatial_object <- Load10X_Spatial(
    data.dir = config$spatial_data_dir,
    bin.size = config$bin_size,
    image.name = config$image_name
  )

  if (!config$assay %in% Assays(spatial_object)) {
    stop(sprintf("The configured assay '%s' was not found.", config$assay))
  }

  DefaultAssay(spatial_object) <- config$assay
  spatial_object
}

plot_tissue_coordinates <- function(spatial_object, config) {
  tissue_plot <- SpatialDimPlot(
    spatial_object,
    plot_segmentations = TRUE,
    alpha = 0.7,
    stroke = 0.04,
    image.scale = "hires"
  ) +
    ggtitle("Spatial tissue overview") +
    NoLegend()

  coordinate_theme <- theme(
    axis.text.x = element_text(size = 8, angle = 45),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black", linetype = "dashed")
  )

  coordinate_plot <- tissue_plot +
    scale_x_continuous(n.breaks = 20) +
    scale_y_reverse(n.breaks = 20) +
    coordinate_theme

  save_plot(coordinate_plot, "tissue_coordinates.pdf", config)
  invisible(coordinate_plot)
}

plot_zoom_region <- function(spatial_object, config) {
  if (is.null(config$zoom_bounds)) {
    message("Zoom plotting skipped because zoom_bounds is not configured.")
    return(invisible(NULL))
  }

  required_bounds <- c("x_min", "x_max", "y_min", "y_max")
  if (!all(required_bounds %in% names(config$zoom_bounds))) {
    stop("zoom_bounds must define x_min, x_max, y_min, and y_max.")
  }

  coordinates <- GetTissueCoordinates(spatial_object, scale = "hires")
  if (!all(c("x", "y") %in% colnames(coordinates))) {
    stop("Tissue coordinates must include x and y columns.")
  }

  coordinates$cell_id <- if ("cell" %in% colnames(coordinates)) {
    coordinates$cell
  } else {
    rownames(coordinates)
  }

  selected_cells <- coordinates %>%
    filter(
      x >= config$zoom_bounds[["x_min"]],
      x <= config$zoom_bounds[["x_max"]],
      y >= config$zoom_bounds[["y_min"]],
      y <= config$zoom_bounds[["y_max"]]
    )

  if (nrow(selected_cells) == 0L) {
    stop("No cells fall within the configured zoom bounds.")
  }

  zoom_object <- subset(spatial_object, cells = selected_cells$cell_id)

  zoom_plot <- SpatialDimPlot(
    zoom_object,
    plot_segmentations = TRUE,
    stroke = 0.1,
    alpha = 0.6,
    image.scale = "hires"
  ) +
    scale_fill_manual(values = config$segmentation_color) +
    ggtitle("Configured tissue region") +
    NoLegend()

  save_plot(zoom_plot, "tissue_zoom.pdf", config)
  invisible(zoom_plot)
}

plot_spatial_counts <- function(spatial_object, config) {
  count_feature <- paste0("nCount_", config$assay)
  if (!count_feature %in% colnames(spatial_object[[]])) {
    stop(sprintf("The spatial count feature '%s' is missing.", count_feature))
  }

  violin_plot <- VlnPlot(
    spatial_object,
    features = count_feature,
    pt.size = 0
  ) +
    theme(axis.text = element_text(size = 12)) +
    NoLegend()

  spatial_count_plot <- SpatialFeaturePlot(
    spatial_object,
    image.scale = "hires",
    features = count_feature
  ) +
    theme(legend.position = "right")

  save_plot(violin_plot, "spatial_count_violin.pdf", config)
  save_plot(spatial_count_plot, "spatial_count_map.pdf", config)
}

run_unsupervised_analysis <- function(spatial_object, config) {
  spatial_object <- NormalizeData(spatial_object)
  spatial_object <- FindVariableFeatures(spatial_object)
  spatial_object <- ScaleData(spatial_object)
  spatial_object <- RunPCA(
    spatial_object,
    assay = config$assay,
    reduction.name = "pca.spatial"
  )
  spatial_object <- FindNeighbors(
    spatial_object,
    assay = config$assay,
    reduction = "pca.spatial",
    dims = config$clustering_dims
  )
  spatial_object <- FindClusters(
    spatial_object,
    cluster.name = "spatial_cluster",
    resolution = config$clustering_resolution
  )
  spatial_object <- RunUMAP(
    spatial_object,
    reduction = "pca.spatial",
    reduction.name = "umap.spatial",
    return.model = TRUE,
    dims = config$clustering_dims
  )

  cluster_plot <- DimPlot(
    spatial_object,
    reduction = "umap.spatial",
    group.by = "spatial_cluster",
    label = TRUE
  ) +
    ggtitle("Unsupervised spatial cell clustering") +
    theme(legend.position = "bottom")

  save_plot(cluster_plot, "spatial_clusters.pdf", config)
  spatial_object
}

transfer_reference_labels <- function(spatial_object, config) {
  reference_object <- load_single_object(config$reference_object_file)

  missing_labels <- setdiff(config$reference_label_columns, colnames(reference_object[[]]))
  if (length(missing_labels) > 0L) {
    stop(
      sprintf(
        "Reference metadata columns are missing: %s",
        paste(missing_labels, collapse = ", ")
      )
    )
  }

  anchors <- FindTransferAnchors(
    reference = reference_object,
    query = spatial_object,
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    dims = config$reduction_dims
  )

  for (label_column in config$reference_label_columns) {
    predictions <- TransferData(
      anchorset = anchors,
      refdata = reference_object[[label_column, drop = TRUE]],
      prediction.assay = FALSE,
      weight.reduction = spatial_object[["pca.spatial"]],
      dims = config$reduction_dims
    )

    metadata_column <- paste0("predicted_", make.names(label_column))
    spatial_object[[metadata_column]] <- predictions$predicted.id
  }

  spatial_object
}

plot_predicted_labels <- function(spatial_object, config) {
  predicted_columns <- paste0(
    "predicted_",
    make.names(config$reference_label_columns)
  )
  missing_columns <- setdiff(predicted_columns, colnames(spatial_object[[]]))
  if (length(missing_columns) > 0L) {
    stop("Predicted labels are not available for plotting.")
  }

  for (metadata_column in predicted_columns) {
    label_values <- sort(unique(na.omit(spatial_object[[metadata_column]][, 1])))
    if (length(label_values) == 0L) {
      warning(sprintf("No labels are available in '%s'.", metadata_column))
      next
    }

    plot <- SpatialDimPlot(
      spatial_object,
      group.by = metadata_column,
      plot_segmentations = TRUE,
      alpha = 0.1,
      stroke = 0.04,
      image.scale = "hires"
    ) +
      labs(fill = metadata_column)

    file_name <- paste0(make.names(metadata_column), "_spatial_plot.pdf")
    save_plot(plot, file_name, config)
  }
}

plot_configured_features <- function(spatial_object, config) {
  if (length(config$features_to_plot) == 0L) {
    message("Feature plotting skipped because features_to_plot is empty.")
    return(invisible(NULL))
  }

  missing_features <- setdiff(
    config$features_to_plot,
    rownames(spatial_object[[config$assay]])
  )
  if (length(missing_features) > 0L) {
    stop(
      sprintf(
        "Configured features are missing: %s",
        paste(missing_features, collapse = ", ")
      )
    )
  }

  for (feature in config$features_to_plot) {
    plot <- SpatialFeaturePlot(
      spatial_object,
      image.scale = "hires",
      features = feature
    ) +
      ggtitle(paste("Expression:", feature))

    file_name <- paste0(make.names(feature), "_spatial_feature.pdf")
    save_plot(plot, file_name, config)
  }
}

main <- function(config) {
  validate_config(config)
  plan(future::multisession, workers = config$workers)
  options(future.globals.maxSize = config$max_memory_gb * 1024^3)

  spatial_object <- load_spatial_object(config)
  plot_tissue_coordinates(spatial_object, config)
  plot_zoom_region(spatial_object, config)
  plot_spatial_counts(spatial_object, config)

  spatial_object <- run_unsupervised_analysis(spatial_object, config)
  spatial_object <- transfer_reference_labels(spatial_object, config)
  plot_predicted_labels(spatial_object, config)
  plot_configured_features(spatial_object, config)

  output_file <- file.path(
    config$output_dir,
    paste0(config$sample_id, "_spatial_cell_segmentation.RData")
  )
  save(spatial_object, file = output_file)
  invisible(spatial_object)
}

spatial_object <- main(config)
