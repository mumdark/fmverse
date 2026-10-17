# Shared downstream result, plotting, and embedding helpers.

fmverse_check_result <- function(result) {
  if (!inherits(result, "fmverse_result")) {
    stop("`result` must be an fmverse_result returned by fmverse().", call. = FALSE)
  }
  if (is.null(result$seurat) || is.null(result$embeddings)) {
    stop("The result must contain both a Seurat object and cell embeddings.", call. = FALSE)
  }
  if (is.null(result$analyses)) result$analyses <- list()
  if (is.null(result$plots)) result$plots <- list()
  if (is.null(result$history)) result$history <- list()
  result
}

fmverse_detach_plot <- function(plot) {
  if (!inherits(plot, "ggplot")) return(plot)
  detach_mapping <- function(mapping) {
    if (!length(mapping)) return(mapping)
    # ggplot2 >= 4.0 stores mappings as an S7 `ggplot2::mapping`, so rebuilding
    # the object with lapply() drops its class and makes the plot invalid. Update
    # each quosure in place through the mapping's `[[<-` method instead.
    for (name in names(mapping)) {
      value <- mapping[[name]]
      if (rlang::is_quosure(value)) {
        mapping[[name]] <- rlang::quo_set_env(value, baseenv())
      }
    }
    mapping
  }
  detach_formula <- function(value) {
    if (inherits(value, "formula")) {
      environment(value) <- baseenv()
      return(value)
    }
    if (is.list(value)) return(lapply(value, detach_formula))
    value
  }
  plot$plot_env <- baseenv()
  plot$mapping <- detach_mapping(plot$mapping)
  for (i in seq_along(plot$layers)) {
    plot$layers[[i]]$mapping <- detach_mapping(plot$layers[[i]]$mapping)
  }
  if (!is.null(plot$facet$params$facets)) {
    plot$facet$params$facets <- detach_formula(plot$facet$params$facets)
  }
  plot
}

fmverse_save_analysis <- function(result, name, analysis, plots = list(), params = list()) {
  result$analyses[[name]] <- analysis
  for (plot_name in names(plots)) {
    result$plots[[plot_name]] <- fmverse_detach_plot(plots[[plot_name]])
  }
  # A chained result may already contain plots built by earlier versions or
  # user code. Detach those too so their mapping environments cannot retain
  # model objects or analysis frames when the result is serialized.
  result$plots <- lapply(result$plots, fmverse_detach_plot)
  result$history[[length(result$history) + 1L]] <- list(
    task = name,
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
    parameters = params
  )
  result
}

fmverse_get_embedding <- function(result, reduction = "scGPT") {
  seurat <- result$seurat
  available <- SeuratObject::Reductions(seurat)
  if (identical(reduction, "scGPT") && !"scGPT" %in% available &&
      !is.null(result$model_name) && result$model_name %in% available) {
    reduction <- result$model_name
  }
  if (reduction %in% available) {
    x <- SeuratObject::Embeddings(seurat[[reduction]])
  } else if (identical(reduction, "scGPT")) {
    x <- result$embeddings
  } else {
    stop("Reduction `", reduction, "` was not found in the Seurat object.", call. = FALSE)
  }
  if (is.null(rownames(x))) stop("Embedding rows must have cell names.", call. = FALSE)
  x
}

fmverse_token_names <- function(result, ids) {
  ids <- as.integer(ids)
  if (identical(result$model_name, "Geneformer") && !is.null(result$inputs$gene_symbols)) {
    labels <- unname(result$inputs$gene_symbols[as.character(ids)])
    fallback <- names(result$inputs$vocab)[match(ids, unname(result$inputs$vocab))]
    labels[is.na(labels) | !nzchar(labels)] <- fallback[is.na(labels) | !nzchar(labels)]
    return(labels)
  }
  names(result$inputs$vocab)[match(ids, unname(result$inputs$vocab))]
}

fmverse_resolve_reduction <- function(result, reduction = "scGPT") {
  available <- SeuratObject::Reductions(result$seurat)
  if (identical(reduction, "scGPT") && !reduction %in% available &&
      identical(result$model_name, "Geneformer") && "Geneformer" %in% available) {
    return("Geneformer")
  }
  reduction
}

fmverse_geneformer_token_id <- function(result, gene) {
  if (methods::is(result, "fmobj")) result <- .fm_view(result, .fm_id(result))
  vocab <- result$inputs$vocab
  token <- unname(vocab[gene])
  if ((!length(token) || is.na(token)) && !is.null(result$inputs$symbol_to_token)) {
    token <- unname(result$inputs$symbol_to_token[gene])
  }
  if (!length(token) || is.na(token)) NA_integer_ else as.integer(token[[1L]])
}

fmverse_check_cells <- function(x, meta) {
  if (!all(rownames(x) %in% rownames(meta))) {
    stop("Embedding cell names do not match Seurat metadata.", call. = FALSE)
  }
  meta[rownames(x), , drop = FALSE]
}

fmverse_with_seed <- function(seed, code) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (!identical(RNGkind(), old_kind)) do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

fmverse_make_umap <- function(seurat, reduction, dims, reduction_name,
                              n_neighbors = 30L, min_dist = 0.3, seed = 1L) {
  reductions <- SeuratObject::Reductions(seurat)
  if (!reduction_name %in% reductions) {
    if (!reduction %in% reductions) {
      stop("Reduction `", reduction, "` is unavailable for UMAP.", call. = FALSE)
    }
    available_dims <- ncol(SeuratObject::Embeddings(seurat[[reduction]]))
    dims <- unique(as.integer(dims))
    dims <- dims[dims >= 1L & dims <= available_dims]
    if (length(dims) < 2L) stop("UMAP requires at least two valid dimensions.", call. = FALSE)
    n_cells <- nrow(SeuratObject::Embeddings(seurat[[reduction]]))
    if (n_cells < 3L) {
      stop("UMAP requires at least three cells.", call. = FALSE)
    }
    # uwot requires n_neighbors < number of observations. Keep the standard
    # default for normal datasets while allowing small fixtures/subsets to run.
    n_neighbors <- min(as.integer(n_neighbors), n_cells - 1L)
    seurat <- Seurat::RunUMAP(
      seurat, reduction = reduction, dims = dims, n.neighbors = as.integer(n_neighbors),
      min.dist = min_dist, seed.use = as.integer(seed), reduction.name = reduction_name,
      verbose = FALSE
    )
  }
  list(seurat = seurat, coordinates = SeuratObject::Embeddings(seurat[[reduction_name]]))
}

fmverse_umap_plot <- function(coordinates, metadata, colour_by, title = NULL,
                              point_size = 0.35) {
  coordinates <- as.matrix(coordinates)
  if (ncol(coordinates) < 2L) stop("UMAP reduction has fewer than two dimensions.")
  metadata <- metadata[rownames(coordinates), , drop = FALSE]
  if (!colour_by %in% colnames(metadata)) stop("Metadata column not found: ", colour_by)
  data <- data.frame(
    cell = rownames(coordinates),
    UMAP_1 = coordinates[, 1L], UMAP_2 = coordinates[, 2L],
    colour = as.factor(metadata[[colour_by]]), stringsAsFactors = FALSE
  )
  ggplot2::ggplot(data, ggplot2::aes(x = UMAP_1, y = UMAP_2, colour = colour)) +
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    ggplot2::labs(title = title, colour = colour_by, x = "UMAP_1", y = "UMAP_2") +
    ggplot2::coord_equal() + ggplot2::theme_classic()
}

fmverse_make_umap_plot <- function(seurat, reduction, dims, reduction_name,
                                   colour_by, title = NULL, seed = 1L) {
  made <- fmverse_make_umap(
    seurat, reduction, dims, reduction_name, seed = seed
  )
  meta <- SeuratObject::FetchData(made$seurat, vars = colour_by)
  names(meta) <- colour_by
  rownames(meta) <- rownames(made$coordinates)
  list(
    seurat = made$seurat,
    plot = fmverse_umap_plot(made$coordinates, meta, colour_by, title),
    coordinates = made$coordinates
  )
}

#' Cluster cells using the fmverse embedding
#'
#' @param result An `fmverse_result`.
#' @param method `"leiden"` or `"louvain"`.
#' @param reduction Seurat reduction to cluster.
#' @param dims Embedding dimensions used to construct the neighbor graph.
#' @param resolution Clustering resolution.
#' @param neighbors Number of nearest neighbors.
#' @param cluster_name Metadata column to write.
#' @param seed Random seed used by UMAP.
#' @return Updated `fmverse_result` with cluster labels and a ggplot UMAP.
#' @export
.fm_cluster_legacy <- function(result, method = c("leiden", "louvain"),
                       reduction = "scGPT", dims = 1:30, resolution = 0.8,
                       neighbors = 20L, cluster_name = "fm_cluster", seed = 1L) {
  result <- fmverse_check_result(result)
  method <- match.arg(method)
  reduction <- fmverse_resolve_reduction(result, reduction)
  if (length(resolution) != 1L || is.na(resolution) || resolution <= 0) {
    stop("`resolution` must be a positive number.", call. = FALSE)
  }
  if (method == "leiden" && !requireNamespace("leidenbase", quietly = TRUE)) {
    stop("Leiden clustering requires the optional package `leidenbase`.", call. = FALSE)
  }
  x <- fmverse_get_embedding(result, reduction)
  dims <- unique(as.integer(dims))
  dims <- dims[dims >= 1L & dims <= ncol(x)]
  if (length(dims) < 2L) stop("At least two valid embedding dimensions are required.")
  nn <- paste0("fm_", method, "_nn")
  snn <- paste0("fm_", method, "_snn")
  seurat <- Seurat::FindNeighbors(
    result$seurat, reduction = reduction, dims = dims, k.param = as.integer(neighbors),
    graph.name = c(nn, snn), verbose = FALSE
  )
  previous_columns <- colnames(seurat[[]])
  seurat <- Seurat::FindClusters(
    seurat, graph.name = snn, algorithm = if (method == "leiden") 4L else 1L,
    resolution = resolution, random.seed = as.integer(seed),
    cluster.name = paste0("fm_tmp_", method, "_", format(resolution, trim = TRUE)),
    verbose = FALSE
  )
  new_columns <- setdiff(colnames(seurat[[]]), previous_columns)
  cluster_column <- paste0("fm_tmp_", method, "_", format(resolution, trim = TRUE))
  if (!cluster_column %in% colnames(seurat[[]])) cluster_column <- tail(new_columns, 1L)
  if (!length(cluster_column)) stop("Seurat did not return cluster labels.")
  cluster_metadata <- data.frame(
    cluster = as.factor(seurat[[cluster_column, drop = TRUE]]),
    row.names = colnames(seurat), stringsAsFactors = FALSE
  )
  names(cluster_metadata) <- cluster_name
  seurat <- SeuratObject::AddMetaData(seurat, metadata = cluster_metadata)
  umap <- fmverse_make_umap_plot(
    seurat, reduction, dims, "umap.fm_cluster", cluster_name,
    paste(result$model_name, method, "clusters"), seed
  )
  result$seurat <- umap$seurat
  result <- fmverse_save_analysis(result, "cluster", list(
    method = method, resolution = resolution, reduction = reduction,
    dims = dims, cluster_column = cluster_name,
    table = table(umap$seurat[[cluster_name, drop = TRUE]])
  ), list(cluster_umap = umap$plot), list(
    method = method, reduction = reduction, dims = dims, resolution = resolution,
    neighbors = neighbors, cluster_name = cluster_name, seed = seed
  ))
  result
}

fmverse_pairwise_distance <- function(query, reference, block_size = 128L) {
  query <- as.matrix(query)
  reference <- as.matrix(reference)
  qnorm <- rowSums(query * query)
  rnorm <- rowSums(reference * reference)
  output <- matrix(NA_real_, nrow(query), nrow(reference))
  for (first in seq.int(1L, nrow(query), by = as.integer(block_size))) {
    ix <- first:min(first + as.integer(block_size) - 1L, nrow(query))
    d2 <- outer(qnorm[ix], rnorm, "+") - 2 * query[ix, , drop = FALSE] %*% t(reference)
    output[ix, ] <- sqrt(pmax(d2, 0))
  }
  output
}

fmverse_knn_group_fraction <- function(embedding, groups, k = 20L, block_size = 128L) {
  embedding <- as.matrix(embedding)
  groups <- as.character(groups)
  if (nrow(embedding) != length(groups) || nrow(embedding) < 2L) return(NA_real_)
  k <- max(1L, min(as.integer(k), nrow(embedding) - 1L))
  squared_norm <- rowSums(embedding * embedding)
  same_fraction <- numeric(nrow(embedding))
  for (first in seq.int(1L, nrow(embedding), by = as.integer(block_size))) {
    index <- first:min(first + as.integer(block_size) - 1L, nrow(embedding))
    distance_squared <- outer(squared_norm[index], squared_norm, "+") -
      2 * embedding[index, , drop = FALSE] %*% t(embedding)
    for (row in seq_along(index)) {
      distance_squared[row, index[[row]]] <- Inf
      nearest <- order(distance_squared[row, ], method = "radix")[seq_len(k)]
      same_fraction[index[[row]]] <- mean(groups[nearest] == groups[index[[row]]], na.rm = TRUE)
    }
  }
  mean(same_fraction, na.rm = TRUE)
}

#' Transfer cell labels from a reference fmverse result
#'
#' @param result Query `fmverse_result`.
#' @param reference Reference `fmverse_result`.
#' @param label Reference metadata column containing labels.
#' @param reduction Shared embedding reduction.
#' @param k Number of reference neighbors.
#' @param predicted_name Query metadata output column.
#' @param truth Optional query metadata column for a confusion plot.
#' @param seed UMAP seed.
#' @return Updated `fmverse_result` with transferred labels and ggplot objects.
#' @export
.fm_annotate_legacy <- function(result, reference, label, reduction = "scGPT", k = 20L,
                        predicted_name = "fm_annotation", truth = NULL, seed = 1L) {
  result <- fmverse_check_result(result)
  reference <- fmverse_check_result(reference)
  if (!identical(result$model_name, reference$model_name)) {
    stop("Query and reference must use the same foundation model.")
  }
  reduction <- fmverse_resolve_reduction(result, reduction)
  reference_reduction <- fmverse_resolve_reduction(reference, reduction)
  query_x <- fmverse_get_embedding(result, reduction)
  reference_x <- fmverse_get_embedding(reference, reference_reduction)
  if (ncol(query_x) != ncol(reference_x)) stop("Query and reference embedding dimensions differ.")
  q_checkpoint <- result$metadata$checkpoint_name
  r_checkpoint <- reference$metadata$checkpoint_name
  if (identical(result$model_name, "Geneformer") && !is.null(q_checkpoint) &&
      !is.null(r_checkpoint) && !identical(q_checkpoint, r_checkpoint)) {
    stop("Query and reference embeddings must use the same Geneformer checkpoint.")
  }
  if (identical(result$model_name, "scGPT") && !is.null(q_checkpoint) && !is.null(r_checkpoint)) {
    q_identity <- fmverse_normalize_scgpt_checkpoint_name(q_checkpoint)
    r_identity <- fmverse_normalize_scgpt_checkpoint_name(r_checkpoint)
    if (!is.na(q_identity) && !is.na(r_identity) && !identical(q_identity, r_identity)) {
      stop("Query and reference embeddings must use the same scGPT checkpoint.")
    }
  }
  if (is.null(rownames(reference$seurat[[]])) || !label %in% colnames(reference$seurat[[]])) {
    stop("Reference label column not found: ", label)
  }
  ref_meta <- reference$seurat[[]][rownames(reference_x), , drop = FALSE]
  ref_labels <- as.character(ref_meta[[label]])
  keep <- !is.na(ref_labels) & nzchar(ref_labels)
  if (!any(keep)) stop("Reference label column contains no usable labels.")
  distance <- fmverse_pairwise_distance(query_x, reference_x[keep, , drop = FALSE])
  k <- max(1L, min(as.integer(k), ncol(distance)))
  ref_labels <- ref_labels[keep]
  predicted <- confidence <- rep(NA, nrow(query_x))
  label_levels <- unique(ref_labels)
  for (i in seq_len(nrow(distance))) {
    nearest <- order(distance[i, ], method = "radix")[seq_len(k)]
    weights <- 1 / pmax(distance[i, nearest], 1e-8)
    vote <- tapply(weights, ref_labels[nearest], sum)
    predicted[[i]] <- names(vote)[which.max(vote)]
    confidence[[i]] <- max(vote) / sum(vote)
  }
  names(predicted) <- names(confidence) <- rownames(query_x)
  predicted_metadata <- data.frame(
    predicted_label = factor(predicted, levels = label_levels),
    predicted_confidence = as.numeric(confidence), row.names = rownames(query_x)
  )
  names(predicted_metadata) <- c(predicted_name, paste0(predicted_name, "_confidence"))
  seurat <- SeuratObject::AddMetaData(result$seurat, predicted_metadata)
  umap <- fmverse_make_umap_plot(
    seurat, reduction, seq_len(min(30L, ncol(query_x))), "umap.fm_annotation",
    predicted_name, paste(result$model_name, "label transfer"), seed
  )
  plots <- list(annotation_umap = umap$plot)
  confusion <- NULL
  if (!is.null(truth)) {
    if (!truth %in% colnames(umap$seurat[[]])) stop("Truth column not found: ", truth)
    confusion <- as.data.frame(table(
      truth = umap$seurat[[truth, drop = TRUE]],
      predicted = umap$seurat[[predicted_name, drop = TRUE]]
    ), stringsAsFactors = FALSE)
    plots$annotation_confusion <- ggplot2::ggplot(
      confusion, ggplot2::aes(x = predicted, y = truth, fill = Freq)
    ) + ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c() +
      ggplot2::labs(x = "Predicted", y = "Reference label", fill = "Cells") +
      ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  result$seurat <- umap$seurat
  result <- fmverse_save_analysis(result, "annotation", list(
    label = label, prediction_column = predicted_name,
    confidence_column = paste0(predicted_name, "_confidence"), confusion = confusion,
    predictions = data.frame(cell = names(predicted), label = unname(predicted),
                             confidence = unname(confidence))
  ), plots, list(label = label, reduction = reduction, k = k,
                 predicted_name = predicted_name, truth = truth, seed = seed))
  result
}

#' Integrate batches using Harmony
#'
#' @param result `fmverse_result` with a batch column in cell metadata.
#' @param batch Metadata column to correct.
#' @param method Currently `"harmony"`.
#' @param reduction Input reduction.
#' @param dims Dimensions used by Harmony and UMAP.
#' @param reduction_name Name of the integrated reduction.
#' @param group_by Optional cell-type column used for the second diagnostic plot.
#' @param seed UMAP seed.
#' @return Updated `fmverse_result` with integrated reduction and ggplots.
#' @export
.fm_integrate_legacy <- function(result, batch, method = "harmony", reduction = "scGPT",
                         dims = 1:30, reduction_name = "fm_harmony",
                         group_by = NULL, seed = 1L) {
  result <- fmverse_check_result(result)
  if (!identical(tolower(method), "harmony")) stop("Only method = 'harmony' is currently supported.")
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("Harmony integration requires the optional package `harmony`.", call. = FALSE)
  }
  reduction <- fmverse_resolve_reduction(result, reduction)
  meta <- result$seurat[[]]
  if (!batch %in% colnames(meta)) stop("Batch metadata column not found: ", batch)
  if (!is.null(group_by) && !group_by %in% colnames(meta)) stop("Metadata column not found: ", group_by)
  original <- fmverse_get_embedding(result, reduction)
  dims <- unique(as.integer(dims))
  dims <- dims[dims >= 1L & dims <= ncol(original)]
  if (length(dims) < 2L) stop("At least two valid dimensions are required for integration.")
  seurat <- harmony::RunHarmony(
    result$seurat, group.by.vars = batch, reduction.use = reduction,
    dims.use = dims, reduction.save = reduction_name, project.dim = FALSE, verbose = FALSE
  )
  batch_labels <- as.character(meta[rownames(original), batch])
  integrated_embedding <- SeuratObject::Embeddings(seurat[[reduction_name]])
  metrics <- list(
    batch_same_neighbor_fraction_before = fmverse_knn_group_fraction(
      original[, dims, drop = FALSE], batch_labels
    ),
    batch_same_neighbor_fraction_after = fmverse_knn_group_fraction(
      integrated_embedding, batch_labels
    )
  )
  if (!is.null(group_by)) {
    cell_labels <- as.character(meta[rownames(original), group_by])
    metrics$cell_type_same_neighbor_fraction_before <- fmverse_knn_group_fraction(
      original[, dims, drop = FALSE], cell_labels
    )
    metrics$cell_type_same_neighbor_fraction_after <- fmverse_knn_group_fraction(
      integrated_embedding, cell_labels
    )
  }
  umap_before <- fmverse_make_umap_plot(
    seurat, reduction, dims, "umap.fm_before_integration", batch,
    "Before Harmony: batch", seed
  )
  seurat <- umap_before$seurat
  umap_after <- fmverse_make_umap_plot(
    seurat, reduction_name, seq_len(min(length(dims), ncol(SeuratObject::Embeddings(seurat[[reduction_name]])))),
    "umap.fm_after_integration", batch, "After Harmony: batch", seed
  )
  seurat <- umap_after$seurat
  plots <- list(integration_batch_before = umap_before$plot,
                integration_batch_after = umap_after$plot)
  if (!is.null(group_by)) {
    before <- fmverse_make_umap_plot(
      seurat, reduction, dims, "umap.fm_before_integration", group_by,
      paste("Before Harmony:", group_by), seed
    )
    after <- fmverse_make_umap_plot(
      seurat, reduction_name, seq_len(min(length(dims), ncol(SeuratObject::Embeddings(seurat[[reduction_name]])))),
      "umap.fm_after_integration", group_by, paste("After Harmony:", group_by), seed
    )
    seurat <- after$seurat
    plots$integration_celltype_before <- before$plot
    plots$integration_celltype_after <- after$plot
  }
  result$seurat <- seurat
  result$analyses$integration <- list(
    method = "harmony", batch = batch, input_reduction = reduction,
    reduction = reduction_name, dims = dims,
    embeddings = integrated_embedding, metrics = metrics
  )
  result <- fmverse_save_analysis(result, "integration", result$analyses$integration,
    plots, list(method = method, batch = batch, reduction = reduction,
                dims = dims, reduction_name = reduction_name, group_by = group_by, seed = seed))
  result
}

#' Find cluster marker genes and return editable Seurat-style plots
#'
#' @param result `fmverse_result` with grouping metadata.
#' @param group_by Metadata column defining groups.
#' @param assay Assay used for differential expression.
#' @param min_pct Minimum expressing-cell fraction.
#' @param logfc_threshold Minimum log fold-change.
#' @param top_n Number of marker genes shown in the dot plot.
#' @return Updated `fmverse_result` with marker table, dot plot and volcano plot.
#' @export
.fm_markers_legacy <- function(result, group_by = "fm_cluster", assay = NULL,
                       min_pct = 0.1, logfc_threshold = 0.25, top_n = 20L) {
  result <- fmverse_check_result(result)
  if (!group_by %in% colnames(result$seurat[[]])) stop("Grouping metadata column not found: ", group_by)
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(result$seurat)
  seurat <- result$seurat
  if (!"data" %in% SeuratObject::Layers(seurat[[assay]])) {
    seurat <- Seurat::NormalizeData(seurat, assay = assay, verbose = FALSE)
  }
  markers <- Seurat::FindAllMarkers(
    seurat, group.by = group_by, assay = assay, only.pos = FALSE,
    min.pct = min_pct, logfc.threshold = logfc_threshold, verbose = FALSE
  )
  if (!nrow(markers)) stop("No markers passed the selected thresholds.", call. = FALSE)
  fc_column <- intersect(c("avg_log2FC", "avg_logFC"), colnames(markers))[[1L]]
  p_column <- if ("p_val_adj" %in% colnames(markers)) "p_val_adj" else "p_val"
  markers$gene <- rownames(markers)
  markers$neg_log10_p <- -log10(pmax(markers[[p_column]], .Machine$double.xmin))
  markers$direction <- ifelse(markers[[fc_column]] >= 0, "up", "down")
  markers$significant <- markers[[p_column]] < 0.05
  volcano <- ggplot2::ggplot(markers, ggplot2::aes(x = .data[[fc_column]], y = neg_log10_p,
                                                  colour = significant)) +
    ggplot2::geom_point(alpha = 0.6, size = 0.8) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey65") +
    ggplot2::geom_vline(xintercept = c(-logfc_threshold, logfc_threshold),
                        linetype = 2, colour = "grey75") +
    # FindAllMarkers names this grouping column `cluster` regardless of the
    # source metadata column supplied through group.by.
    ggplot2::facet_wrap(~cluster, scales = "free") +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
    ggplot2::labs(x = "Average log fold change", y = paste0("-log10(", p_column, ")"),
                  colour = "p < 0.05", title = "Differential expression") + ggplot2::theme_classic()
  top_genes <- unique(utils::head(markers$gene[order(markers[[p_column]], -abs(markers[[fc_column]]))],
                                 as.integer(top_n)))
  top_genes <- intersect(top_genes, rownames(SeuratObject::LayerData(
    seurat, assay = assay, layer = "data"
  )))
  expression <- SeuratObject::LayerData(seurat, assay = assay, layer = "data")
  groups <- as.character(seurat[[group_by, drop = TRUE]])
  dot_data <- do.call(rbind, lapply(unique(groups), function(group) {
    cells <- which(groups == group)
    do.call(rbind, lapply(top_genes, function(gene) {
      values <- as.numeric(expression[gene, cells, drop = TRUE])
      data.frame(group = group, gene = gene, average = mean(values),
                 percent = mean(values > 0), stringsAsFactors = FALSE)
    }))
  }))
  dot <- ggplot2::ggplot(dot_data, ggplot2::aes(x = gene, y = group,
                                                size = percent, colour = average)) +
    ggplot2::geom_point() + ggplot2::scale_colour_viridis_c() +
    ggplot2::labs(title = "Top differential-expression markers", x = NULL, y = group_by,
                  colour = "Average expression", size = "Percent expressed") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  result$seurat <- seurat
  result <- fmverse_save_analysis(result, "markers", list(
    table = markers, group_by = group_by, assay = assay, top_genes = top_genes
  ), list(markers_volcano = volcano, markers_dotplot = dot),
  list(group_by = group_by, assay = assay, min_pct = min_pct,
       logfc_threshold = logfc_threshold, top_n = top_n))
  result
}

#' Run Gene Ontology enrichment on differential-expression markers
#'
#' @param result `fmverse_result` with `fm_markers()` results.
#' @param markers Optional marker table; defaults to the stored table.
#' @param ontology GO ontology (`BP`, `MF`, or `CC`).
#' @param pvalue_cutoff Enrichment p-value cutoff.
#' @param top_n Maximum number of pathways per group shown in the plots.
#' @return Updated `fmverse_result` with enrichment table and ggplots.
#' @export
.fm_enrich_legacy <- function(result, markers = NULL, ontology = "BP", pvalue_cutoff = 0.05,
                      top_n = 5L) {
  result <- fmverse_check_result(result)
  top_n <- as.integer(top_n)
  if (length(top_n) != 1L || is.na(top_n) || top_n < 1L) {
    stop("`top_n` must be a positive integer.", call. = FALSE)
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("GO enrichment requires optional packages `clusterProfiler` and `org.Hs.eg.db`.", call. = FALSE)
  }
  if (is.null(markers)) markers <- result$analyses$markers$table
  if (is.null(markers) || !nrow(markers)) stop("Provide markers or run fm_markers() first.")
  if (!"gene" %in% colnames(markers)) markers$gene <- rownames(markers)
  p_col <- if ("p_val_adj" %in% colnames(markers)) "p_val_adj" else "p_val"
  groups <- if ("cluster" %in% colnames(markers)) unique(as.character(markers$cluster)) else "all"
  all_enrichment <- list()
  for (group in groups) {
    selected <- if (group == "all") markers else markers[as.character(markers$cluster) == group, , drop = FALSE]
    genes <- unique(selected$gene[!is.na(selected[[p_col]]) & selected[[p_col]] < pvalue_cutoff])
    if (length(genes) < 2L) next
    enriched <- clusterProfiler::enrichGO(
      gene = genes, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL",
      ont = ontology, pAdjustMethod = "BH", pvalueCutoff = pvalue_cutoff,
      readable = TRUE
    )
    tab <- as.data.frame(enriched)
    if (nrow(tab)) {
      tab$cluster <- group
      all_enrichment[[group]] <- tab
    }
  }
  table <- if (length(all_enrichment)) do.call(rbind, all_enrichment) else data.frame()
  if (!nrow(table)) {
    empty_data <- data.frame(x = 0, y = 0, label = "No enriched GO terms at this cutoff")
    empty_plot <- ggplot2::ggplot(empty_data, ggplot2::aes(x, y, label = label)) +
      ggplot2::geom_text() + ggplot2::labs(title = paste("GO", ontology, "enrichment")) +
      ggplot2::theme_void()
    return(fmverse_save_analysis(result, "enrichment", list(
      table = table, ontology = ontology,
      note = "No terms passed the requested p-value cutoff."
    ), list(enrichment_dotplot = empty_plot, enrichment_barplot = empty_plot),
    list(ontology = ontology, pvalue_cutoff = pvalue_cutoff, top_n = top_n)))
  }
  # Retain all terms in the analysis table, but bound the number displayed in
  # fixed-size plots so long GO labels do not overlap.
  plot_n <- min(top_n, max(1L, floor(30L / length(unique(table$cluster)))))
  plot_table <- do.call(rbind, lapply(split(table, table$cluster), function(group_table) {
    group_table <- group_table[order(group_table$p.adjust, group_table$pvalue), , drop = FALSE]
    utils::head(group_table, plot_n)
  }))
  wrap_terms <- function(x) vapply(as.character(x), function(term) {
    paste(strwrap(term, width = 46L), collapse = "\n")
  }, character(1))
  plot_table$Description <- wrap_terms(plot_table$Description)
  plot_table$Description <- factor(plot_table$Description,
    levels = rev(unique(plot_table$Description[order(plot_table$p.adjust)])))
  dot <- ggplot2::ggplot(plot_table, ggplot2::aes(x = cluster, y = Description,
                                            size = Count, colour = p.adjust)) +
    ggplot2::geom_point(alpha = 0.9) + ggplot2::scale_colour_viridis_c(direction = -1) +
    ggplot2::labs(x = "Group", y = NULL, colour = "Adjusted p", size = "Genes",
                  title = paste("GO", ontology, "enrichment"),
                  subtitle = paste("Top", plot_n, "terms per group")) + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
  bar_n <- min(3L, plot_n)
  bar_table <- do.call(rbind, lapply(split(table, table$cluster), function(group_table) {
    group_table <- group_table[order(group_table$p.adjust, group_table$pvalue), , drop = FALSE]
    utils::head(group_table, bar_n)
  }))
  bar_table$Description <- wrap_terms(bar_table$Description)
  bar_table$pathway <- paste(bar_table$Description, paste0("(", bar_table$cluster, ")"))
  bar_table$pathway <- factor(bar_table$pathway,
    levels = unique(bar_table$pathway[order(-log10(bar_table$p.adjust))]))
  bar <- ggplot2::ggplot(bar_table,
                         ggplot2::aes(x = pathway, y = -log10(p.adjust), fill = cluster)) +
    ggplot2::geom_col(position = "dodge") + ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "-log10(adjusted p)", title = "Top enriched pathways",
                  subtitle = paste("Up to", bar_n, "terms per group"), fill = "Group") +
    ggplot2::theme_classic() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
  fmverse_save_analysis(result, "enrichment", list(table = table, ontology = ontology),
                        list(enrichment_dotplot = dot, enrichment_barplot = bar),
                        list(ontology = ontology, pvalue_cutoff = pvalue_cutoff, top_n = top_n))
}

#' Infer pseudotime using Slingshot
#'
#' @param result `fmverse_result` with a cluster column.
#' @param group_by Cluster metadata column.
#' @param root Starting cluster.
#' @param reduction Embedding reduction.
#' @param dims Dimensions passed to Slingshot.
#' @param seed UMAP seed for the visualization.
#' @return Updated `fmverse_result` with pseudotime and plots.
#' @export
.fm_trajectory_legacy <- function(result, group_by = "fm_cluster", root,
                          reduction = "scGPT", dims = 1:20, seed = 1L) {
  result <- fmverse_check_result(result)
  reduction <- fmverse_resolve_reduction(result, reduction)
  if (!requireNamespace("slingshot", quietly = TRUE)) {
    stop("Trajectory analysis requires the optional package `slingshot`.", call. = FALSE)
  }
  if (missing(root) || length(root) != 1L || is.na(root)) stop("Provide one root cluster.")
  if (!group_by %in% colnames(result$seurat[[]])) stop("Cluster metadata column not found: ", group_by)
  x <- fmverse_get_embedding(result, reduction)
  dims <- unique(as.integer(dims))
  dims <- dims[dims >= 1L & dims <= ncol(x)]
  if (length(dims) < 2L) stop("At least two valid dimensions are required.")
  meta <- result$seurat[[]][rownames(x), , drop = FALSE]
  clusters <- factor(meta[[group_by]])
  if (!as.character(root) %in% levels(clusters)) stop("Root cluster is not present in `group_by`.")
  sds <- slingshot::slingshot(x[, dims, drop = FALSE], clusterLabels = clusters,
                              start.clus = as.character(root))
  pseudotime <- slingshot::slingPseudotime(sds)
  if (is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime, ncol = 1L)
  rownames(pseudotime) <- rownames(x)
  colnames(pseudotime) <- paste0("fm_pseudotime_", seq_len(ncol(pseudotime)))
  seurat <- SeuratObject::AddMetaData(result$seurat, as.data.frame(pseudotime))
  umap <- fmverse_make_umap(seurat, reduction, dims, "umap.fm_trajectory", seed = seed)
  seurat <- umap$seurat
  pseudo_long <- data.frame(cell = rownames(pseudotime),
                            pseudotime = rowMeans(pseudotime, na.rm = TRUE),
                            stringsAsFactors = FALSE)
  pseudo_long <- pseudo_long[is.finite(pseudo_long$pseudotime), , drop = FALSE]
  coords <- umap$coordinates[pseudo_long$cell, , drop = FALSE]
  pseudo_long$UMAP_1 <- coords[, 1L]
  pseudo_long$UMAP_2 <- coords[, 2L]
  ptime_plot <- ggplot2::ggplot(pseudo_long, ggplot2::aes(UMAP_1, UMAP_2, colour = pseudotime)) +
    ggplot2::geom_point(size = 0.35) + ggplot2::scale_colour_viridis_c() +
    ggplot2::labs(title = "Slingshot pseudotime", colour = "Pseudotime",
                  x = "UMAP_1", y = "UMAP_2") +
    ggplot2::coord_equal() + ggplot2::theme_classic()
  display_coordinates <- umap$coordinates[rownames(x), , drop = FALSE]
  # Slingshot's fitted curves live in the analysis space. Fit a second curve
  # in the 2D display space for a faithful overlay; projecting high-dimensional
  # curve coordinates to UMAP with nearest cells can create misleading jumps.
  display_sds <- slingshot::slingshot(
    display_coordinates[, 1:2, drop = FALSE], clusterLabels = clusters,
    start.clus = as.character(root)
  )
  curves <- as.data.frame(sl_curve_points(display_sds))
  traj_data <- data.frame(x = display_coordinates[, 1L], y = display_coordinates[, 2L],
                          cluster = clusters)
  show_cluster_legend <- nlevels(clusters) <= 10L
  traj_plot <- if (show_cluster_legend) {
    ggplot2::ggplot(traj_data, ggplot2::aes(x = x, y = y, colour = cluster)) +
      ggplot2::geom_point(size = 0.35, alpha = 0.45) +
      ggplot2::labs(colour = group_by)
  } else {
    ggplot2::ggplot(traj_data, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(size = 0.35, alpha = 0.45, colour = "grey55") +
      ggplot2::labs(subtitle = paste("Cells shown in grey;", nlevels(clusters),
                                     "groups omitted from legend"))
  }
  traj_plot <- traj_plot +
    ggplot2::geom_path(data = curves, ggplot2::aes(x = x, y = y, group = lineage),
                       inherit.aes = FALSE, linewidth = 0.65, alpha = 0.7,
                       colour = "grey20") +
    ggplot2::labs(title = "Slingshot lineage curves on UMAP",
                  x = "UMAP_1", y = "UMAP_2") +
    ggplot2::coord_equal() + ggplot2::theme_classic()
  result$seurat <- seurat
  fmverse_save_analysis(result, "trajectory", list(
    root = as.character(root), group_by = group_by, pseudotime = pseudotime,
    lineages = slingshot::slingLineages(sds)
  ), list(trajectory_pseudotime = ptime_plot, trajectory_curves = traj_plot),
  list(group_by = group_by, root = root, reduction = reduction, dims = dims, seed = seed))
}

sl_curve_points <- function(sds) {
  curves <- slingshot::slingCurves(sds)
  pieces <- lapply(seq_along(curves), function(i) {
    curve_points <- as.matrix(curves[[i]]$s)
    curve_points <- curve_points[stats::complete.cases(curve_points), , drop = FALSE]
    if (!nrow(curve_points) || ncol(curve_points) < 2L) return(NULL)
    data.frame(x = curve_points[, 1L], y = curve_points[, 2L],
               lineage = paste0("lineage_", i))
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(data.frame(x = numeric(), y = numeric(), lineage = character()))
  do.call(rbind, pieces)
}

fmverse_ensure_model <- function(result) {
  # Use exact list indexing here: a saved result may omit `model`, and `$model`
  # would partially match `model_name` ("scGPT") and mistake it for a module.
  if ("model" %in% names(result) && !is.null(result[["model"]])) {
    model <- result[["model"]]
    model$to(device = torch::torch_device("cpu"))
    model$eval()
    return(model)
  }
  checkpoint <- result$metadata$checkpoint
  if (is.null(checkpoint) || !file.exists(checkpoint)) {
    stop("This result has no live model; reload it with fmverse() or provide its original checkpoint.")
  }
  if (identical(result$model_name, "Geneformer")) {
    config <- result$metadata$config
    if (is.null(config)) {
      config_path <- file.path(dirname(checkpoint), "config.json")
      if (!file.exists(config_path)) stop("Geneformer config.json is missing beside checkpoint.")
      config <- jsonlite::fromJSON(config_path, simplifyVector = TRUE)
    }
    model <- fmverse_geneformer_build_model(config)
    result[["model"]] <- fmverse_load_geneformer_weights(model, checkpoint)$model
    result[["model"]]$to(device = torch::torch_device("cpu"))
    result[["model"]]$eval()
    return(result[["model"]])
  }
  args_file <- file.path(dirname(checkpoint), "args.json")
  vocab_file <- file.path(dirname(checkpoint), "vocab.json")
  if (!file.exists(args_file) || !file.exists(vocab_file)) stop("Checkpoint args.json/vocab.json is missing.")
  config <- jsonlite::fromJSON(args_file, simplifyVector = TRUE)
  vocab <- jsonlite::fromJSON(vocab_file, simplifyVector = TRUE)
  model <- fmverse_scgpt_model(
    ntoken = length(vocab), pad_id = unname(vocab[["<pad>"]]),
    d_model = config$embsize, nhead = config$nheads, d_hid = config$d_hid,
    nlayers = config$nlayers, dropout = config$dropout
  )
  result[["model"]] <- fmverse_load_scgpt_weights(model, checkpoint,
    nlayers_cls = config$n_layers_cls)$model
  result[["model"]]$to(device = torch::torch_device("cpu"))
  result[["model"]]$eval()
  result[["model"]]
}

fmverse_forward_inputs <- function(result, indices, return_expression = FALSE,
                                   return_attention = FALSE, return_layers = FALSE,
                                   max_length = NULL) {
  if (methods::is(result, "fmobj")) result <- .fm_view(result, .fm_id(result))
  if (is.null(result$inputs)) {
    stop("Token inputs were not stored in this result. Re-run fmverse() with the current package version.")
  }
  model <- fmverse_ensure_model(result)
  model$to(device = torch::torch_device("cpu"))
  model$eval()
  inputs <- result$inputs
  indices <- as.integer(indices)
  # `which` is flattened; use the last occupied position in any selected row.
  active <- if (identical(result$model_name, "Geneformer")) {
    inputs$attention_mask[indices, , drop = FALSE]
  } else {
    !inputs$pad_mask[indices, , drop = FALSE]
  }
  seq_length <- max(rowSums(active), 1L)
  if (!is.null(max_length)) seq_length <- min(seq_length, as.integer(max_length))
  gene_ids <- torch::torch_tensor(inputs$gene_ids[indices, seq_len(seq_length), drop = FALSE],
                                  dtype = torch::torch_long(), device = torch::torch_device("cpu"))
  if (identical(result$model_name, "Geneformer")) {
    attention_mask <- torch::torch_tensor(
      inputs$attention_mask[indices, seq_len(seq_length), drop = FALSE],
      dtype = torch::torch_long(), device = torch::torch_device("cpu")
    )
    output <- torch::with_no_grad(model(
      input_ids = gene_ids, attention_mask = attention_mask,
      return_layers = return_layers, return_attention = return_attention,
      return_logits = return_expression, embedding_layer = -1L
    ))
    return(list(output = output, sequence_length = seq_length, inputs = inputs, model = model,
                attention_mask = attention_mask))
  }
  values <- torch::torch_tensor(inputs$values[indices, seq_len(seq_length), drop = FALSE],
                               dtype = torch::torch_float(), device = torch::torch_device("cpu"))
  pad_mask <- torch::torch_tensor(inputs$pad_mask[indices, seq_len(seq_length), drop = FALSE],
                                 dtype = torch::torch_bool(), device = torch::torch_device("cpu"))
  output <- torch::with_no_grad(model(
    gene_ids, values, pad_mask, return_expression = return_expression,
    return_attention = return_attention
  ))
  list(output = output, sequence_length = seq_length, inputs = inputs, model = model)
}

#' Reconstruct expression values with the scGPT expression decoder
#'
#' @param result `fmverse_result` with checkpoint model and saved token inputs.
#' @param cells Optional cell barcodes to reconstruct.
#' @param max_cells Maximum number of cells for the diagnostic plot.
#' @param seed Seed for deterministic cell subsampling.
#' @return Updated `fmverse_result` with predictions and an editable ggplot.
#' @export
.fm_reconstruct_legacy <- function(result, cells = NULL, max_cells = 128L, seed = 1L) {
  result <- fmverse_check_result(result)
  if (!result$model_name %in% c("scGPT", "Geneformer")) {
    stop("`fm_reconstruct()` currently supports scGPT and Geneformer only; ",
      result$model_name, " has no adapted reconstruction decoder.", call. = FALSE)
  }
  if (identical(result$model_name, "Geneformer")) {
    return(fmverse_geneformer_reconstruct(result, cells, max_cells, seed))
  }
  if (is.null(result$inputs)) stop("Token inputs are absent; re-run fmverse() with the current version.")
  selected <- if (is.null(cells)) result$inputs$cell_names else as.character(cells)
  if (any(!selected %in% result$inputs$cell_names)) stop("Some selected cells are absent from saved tokens.")
  if (length(selected) > max_cells) {
    set.seed(as.integer(seed)); selected <- sample(selected, as.integer(max_cells))
  }
  indices <- match(selected, result$inputs$cell_names)
  forward <- fmverse_forward_inputs(result, indices, return_expression = TRUE)
  result$model <- forward$model
  prediction <- as.matrix(torch::as_array(forward$output$expression_prediction$cpu()))
  observed <- result$inputs$values[indices, seq_len(forward$sequence_length), drop = FALSE]
  ids <- result$inputs$gene_ids[indices, seq_len(forward$sequence_length), drop = FALSE]
  pad <- result$inputs$pad_mask[indices, seq_len(forward$sequence_length), drop = FALSE]
  keep <- !pad & ids != unname(result$inputs$vocab[["<cls>"]])
  token_names <- names(result$inputs$vocab)[match(as.integer(ids), unname(result$inputs$vocab))]
  comparison <- data.frame(
    gene = token_names[keep], observed = observed[keep], predicted = prediction[keep],
    cell = matrix(rep(selected, times = ncol(ids)), nrow = nrow(ids),
                  ncol = ncol(ids))[keep], stringsAsFactors = FALSE
  )
  comparison <- comparison[is.finite(comparison$observed) & is.finite(comparison$predicted), ]
  common_range <- range(c(comparison$observed, comparison$predicted), finite = TRUE)
  plot <- ggplot2::ggplot(comparison, ggplot2::aes(observed, predicted)) +
    ggplot2::geom_point(alpha = 0.16, size = 0.55) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         colour = "firebrick", linewidth = 0.8) +
    ggplot2::labs(x = "Input binned value", y = "Decoder prediction",
                  title = "scGPT expression reconstruction",
                  subtitle = "Dashed line: identity; red line: linear trend") +
    ggplot2::coord_equal(xlim = common_range, ylim = common_range, expand = FALSE) +
    ggplot2::theme_classic()
  result <- fmverse_save_analysis(result, "reconstruct", list(
    comparison = comparison, cells = selected, n_tokens = nrow(comparison),
    interpretation = "Decoder predictions are token-level model outputs, not raw-count estimates."
  ), list(reconstruction_scatter = plot), list(cells = selected, max_cells = max_cells, seed = seed))
  result
}

fmverse_geneformer_reconstruct <- function(result, cells = NULL, max_cells = 128L,
                                            seed = 1L, mask_fraction = 0.15,
                                            max_masked_tokens = 32L) {
  if (is.null(result$inputs$attention_mask) || is.null(result$inputs$vocab[["<mask>"]])) {
    stop("Geneformer masked-token inputs are missing; re-run fmverse() with the current version.")
  }
  selected <- if (is.null(cells)) result$inputs$cell_names else as.character(cells)
  if (any(!selected %in% result$inputs$cell_names)) stop("Some selected cells are absent from saved tokens.")
  if (length(selected) > max_cells) selected <- fmverse_with_seed(seed, sample(selected, max_cells))
  indices <- match(selected, result$inputs$cell_names)
  model <- fmverse_ensure_model(result)
  result$model <- model
  mask_id <- as.integer(result$inputs$vocab[["<mask>"]])
  reconstruction_rows <- vector("list", length(indices))
  for (i in seq_along(indices)) {
    row <- indices[[i]]
    length_i <- sum(result$inputs$attention_mask[row, ])
    if (length_i < 1L) next
    sequence <- result$inputs$gene_ids[row, seq_len(length_i), drop = FALSE]
    eligible <- seq_len(length_i)
    if (model$cls_token_id >= 0L) eligible <- eligible[sequence[1L, eligible] != model$cls_token_id]
    if (model$eos_token_id >= 0L) eligible <- eligible[sequence[1L, eligible] != model$eos_token_id]
    if (!length(eligible)) next
    mask_n <- min(as.integer(max_masked_tokens), max(1L, floor(length(eligible) * mask_fraction)))
    positions <- fmverse_with_seed(seed + i, eligible[sample.int(length(eligible), mask_n)])
    masked_sequence <- sequence
    masked_sequence[1L, positions] <- mask_id
    input_ids <- torch::torch_tensor(masked_sequence, dtype = torch::torch_long(),
                                     device = torch::torch_device("cpu"))
    attention_mask <- torch::torch_ones_like(input_ids)
    output <- torch::with_no_grad(model(
      input_ids = input_ids, attention_mask = attention_mask,
      return_layers = TRUE, embedding_layer = -1L
    ))
    # Only run the large tied MLM projection for held-out positions. Projecting
    # every token in a 4,096-gene sequence wastes most of the checkpoint's
    # vocabulary-matrix compute during this diagnostic.
    selected_hidden <- output$hidden_states[[as.integer(model$num_hidden_layers)]]
    selected_hidden <- selected_hidden[, as.integer(positions), , drop = FALSE]
    logits_tensor <- torch::with_no_grad(model$cls(
      selected_hidden, model$bert$embeddings$word_embeddings$weight
    ))
    logits_array <- torch::as_array(logits_tensor$cpu())
    logits <- matrix(logits_array[1L, , ], nrow = length(positions),
                     ncol = length(result$inputs$vocab))
    true_ids <- as.integer(sequence[1L, positions])
    predicted_ids <- integer(length(positions))
    true_probability <- numeric(length(positions))
    true_rank <- integer(length(positions))
    for (j in seq_along(positions)) {
      scores <- logits[j, ]
      predicted_ids[[j]] <- which.max(scores) - 1L
      shifted <- scores - max(scores)
      probabilities <- exp(shifted) / sum(exp(shifted))
      true_probability[[j]] <- probabilities[true_ids[[j]] + 1L]
      true_rank[[j]] <- sum(scores > scores[true_ids[[j]] + 1L]) + 1L
    }
    reconstruction_rows[[i]] <- data.frame(
      cell = selected[[i]], position = positions, true_token_id = true_ids,
      true_gene = fmverse_token_names(result, true_ids),
      predicted_token_id = predicted_ids,
      predicted_gene = fmverse_token_names(result, predicted_ids),
      true_probability = true_probability, true_rank = true_rank,
      top1_correct = predicted_ids == true_ids, stringsAsFactors = FALSE
    )
    rm(output, selected_hidden, logits_tensor, logits)
  }
  comparison <- do.call(rbind, Filter(Negate(is.null), reconstruction_rows))
  if (is.null(comparison) || !nrow(comparison)) stop("No masked Geneformer tokens were available.")
  summary <- data.frame(
    metric = c("masked_tokens", "top1_accuracy", "median_true_probability", "median_true_rank"),
    value = c(nrow(comparison), mean(comparison$top1_correct),
              stats::median(comparison$true_probability), stats::median(comparison$true_rank))
  )
  plot <- ggplot2::ggplot(comparison, ggplot2::aes(x = position, y = true_probability)) +
    ggplot2::geom_point(alpha = 0.35, size = 0.7, colour = "#2C7FB8") +
    ggplot2::geom_smooth(method = "loess", se = FALSE, colour = "#D95F0E", linewidth = 0.7) +
    ggplot2::labs(
      title = "Geneformer masked-token reconstruction",
      subtitle = sprintf("Top-1 accuracy %.3f; masked genes are held out before prediction",
                         mean(comparison$top1_correct)),
      x = "Rank position in input sequence", y = "Probability assigned to true gene"
    ) + ggplot2::theme_classic()
  fmverse_save_analysis(result, "reconstruct", list(
    comparison = comparison, metrics = summary, cells = selected,
    interpretation = "Masked-token prediction diagnostic for the Geneformer MLM head; not raw-count reconstruction."
  ), list(reconstruction_scatter = plot), list(
    cells = selected, max_cells = max_cells, seed = seed,
    mask_fraction = mask_fraction, max_masked_tokens = max_masked_tokens,
    embedding_layer = -1L
  ))
}

#' Perform in-silico gene perturbation with the scGPT decoder
#'
#' @param result `fmverse_result` with saved token inputs.
#' @param genes Gene symbols to perturb.
#' @param cells Optional cells to include.
#' @param max_cells Maximum cells used per perturbation.
#' @param top_n Number of genes shown per perturbation.
#' @param seed Random seed for subsampling.
#' @return Updated `fmverse_result` with a predicted gene-effect table and heatmap.
#' @export
.fm_perturb_legacy <- function(result, genes, cells = NULL, max_cells = 128L,
                       top_n = 20L, seed = 1L) {
  result <- fmverse_check_result(result)
  if (!result$model_name %in% c("scGPT", "Geneformer")) {
    stop("`fm_perturb()` currently supports scGPT and Geneformer only; ",
      result$model_name, " has no adapted perturbation interface.", call. = FALSE)
  }
  if (identical(result$model_name, "Geneformer")) {
    return(fmverse_geneformer_perturb(result, genes, cells, max_cells, top_n, seed))
  }
  if (is.null(result$inputs)) stop("Token inputs are absent; re-run fmverse() with the current version.")
  genes <- unique(as.character(genes))
  if (!length(genes) || anyNA(genes)) stop("Provide one or more gene symbols.")
  selected <- if (is.null(cells)) result$inputs$cell_names else as.character(cells)
  if (any(!selected %in% result$inputs$cell_names)) stop("Some selected cells are absent from saved tokens.")
  if (length(selected) > max_cells) {
    set.seed(as.integer(seed)); selected <- sample(selected, as.integer(max_cells))
  }
  indices <- match(selected, result$inputs$cell_names)
  result$model <- fmverse_ensure_model(result)
  base <- fmverse_forward_inputs(result, indices, return_expression = TRUE)
  base_prediction <- as.matrix(torch::as_array(base$output$expression_prediction$cpu()))
  ids <- result$inputs$gene_ids[indices, seq_len(base$sequence_length), drop = FALSE]
  original_values <- result$inputs$values[indices, seq_len(base$sequence_length), drop = FALSE]
  pad <- result$inputs$pad_mask[indices, seq_len(base$sequence_length), drop = FALSE]
  token_names <- matrix(
    names(result$inputs$vocab)[match(as.integer(ids), unname(result$inputs$vocab))],
    nrow = nrow(ids), ncol = ncol(ids)
  )
  outputs <- list()
  cells_by_perturbation <- list()
  for (gene in genes) {
    gene_id <- result$inputs$vocab[[gene]]
    if (is.null(gene_id) || is.na(gene_id)) {
      warning("Gene not present in scGPT vocabulary: ", gene, call. = FALSE)
      next
    }
    target_cells <- which(rowSums(ids == as.integer(gene_id)) > 0L)
    if (!length(target_cells)) {
      warning("No selected cell contains the target token `", gene,
              "`; it was skipped for perturbation.", call. = FALSE)
      next
    }
    cells_by_perturbation[[gene]] <- selected[target_cells]
    gene_ids <- ids[target_cells, , drop = FALSE]
    gene_values <- original_values[target_cells, , drop = FALSE]
    gene_pad <- pad[target_cells, , drop = FALSE]
    gene_base_prediction <- base_prediction[target_cells, , drop = FALSE]
    gene_token_names <- token_names[target_cells, , drop = FALSE]
    modified <- gene_values
    modified[gene_ids == as.integer(gene_id)] <- 0
    altered <- torch::torch_tensor(modified, dtype = torch::torch_float(), device = torch::torch_device("cpu"))
    gene_tensor <- torch::torch_tensor(gene_ids, dtype = torch::torch_long(), device = torch::torch_device("cpu"))
    pad_tensor <- torch::torch_tensor(gene_pad, dtype = torch::torch_bool(), device = torch::torch_device("cpu"))
    perturbed_output <- torch::with_no_grad(result$model(
      gene_tensor, altered, pad_tensor, return_expression = TRUE
    ))
    perturbed <- as.matrix(torch::as_array(perturbed_output$expression_prediction$cpu()))
    effect <- perturbed - gene_base_prediction
    target_positions <- gene_ids == as.integer(gene_id)
    effect[target_positions] <- NA_real_
    per_gene <- tapply(as.numeric(effect), as.vector(gene_token_names), mean, na.rm = TRUE)
    per_gene <- per_gene[is.finite(per_gene)]
    per_gene <- per_gene[setdiff(names(per_gene), c("<pad>", "<cls>", gene))]
    if (length(per_gene)) {
      selected_effects <- head(per_gene[order(abs(per_gene), decreasing = TRUE)], as.integer(top_n))
      outputs[[gene]] <- data.frame(perturbation = gene, gene = names(selected_effects),
                                    predicted_effect = as.numeric(selected_effects),
                                    stringsAsFactors = FALSE)
    }
  }
  if (!length(outputs)) stop("No requested genes were found in the scGPT vocabulary or token inputs.")
  effects <- do.call(rbind, outputs)
  effects$gene <- factor(effects$gene, levels = unique(effects$gene[order(effects$predicted_effect)]))
  effect_limit <- max(abs(effects$predicted_effect), na.rm = TRUE)
  if (!is.finite(effect_limit) || effect_limit == 0) effect_limit <- 1
  cell_counts <- vapply(cells_by_perturbation, length, integer(1))
  subtitle <- paste0("Cells carrying perturbed token: ",
                     paste0(names(cell_counts), " n=", cell_counts, collapse = "; "))
  heatmap <- ggplot2::ggplot(effects, ggplot2::aes(x = perturbation, y = gene,
                                                   fill = predicted_effect)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F2F2F2", high = "#B2182B",
                                  midpoint = 0, limits = c(-effect_limit, effect_limit),
                                  name = "Decoder delta") +
    ggplot2::labs(x = "Perturbed gene", y = "Predicted affected gene",
                  title = "In-silico perturbation effects", subtitle = subtitle) +
    ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  fmverse_save_analysis(result, "perturbation", list(
    effects = effects, cells = selected, cells_by_perturbation = cells_by_perturbation,
    interpretation = "Exploratory decoder response after zeroing the perturbed token value; not a causal effect estimate."
  ), list(perturbation_heatmap = heatmap), list(genes = genes, cells = selected,
                                                max_cells = max_cells, top_n = top_n, seed = seed))
}

fmverse_geneformer_perturb <- function(result, genes, cells = NULL, max_cells = 128L,
                                       top_n = 20L, seed = 1L) {
  if (is.null(result$inputs$attention_mask) || is.null(result$inputs$vocab[["<mask>"]])) {
    stop("Geneformer token inputs are missing; re-run fmverse() with the current version.")
  }
  genes <- unique(as.character(genes))
  if (!length(genes) || anyNA(genes)) stop("Provide one or more gene symbols or Ensembl IDs.")
  selected <- if (is.null(cells)) result$inputs$cell_names else as.character(cells)
  if (any(!selected %in% result$inputs$cell_names)) stop("Some selected cells are absent from saved tokens.")
  if (length(selected) > max_cells) selected <- fmverse_with_seed(seed, sample(selected, max_cells))
  selected_indices <- match(selected, result$inputs$cell_names)
  model <- fmverse_ensure_model(result)
  result$model <- model
  mask_id <- as.integer(result$inputs$vocab[["<mask>"]])
  hidden_index <- as.integer(model$num_hidden_layers)
  outputs <- list()
  cells_by_perturbation <- list()
  for (gene in genes) {
    target_id <- fmverse_geneformer_token_id(result, gene)
    if (is.na(target_id)) {
      warning("Gene not present in Geneformer vocabulary: ", gene, call. = FALSE)
      next
    }
    has_target <- vapply(selected_indices, function(i) {
      any(result$inputs$gene_ids[i, result$inputs$attention_mask[i, ]] == target_id)
    }, logical(1))
    target_indices <- selected_indices[has_target]
    if (!length(target_indices)) {
      warning("No selected cell contains target token `", gene, "`; skipped.", call. = FALSE)
      next
    }
    cells_by_perturbation[[gene]] <- result$inputs$cell_names[target_indices]
    gene_effects <- list()
    for (cell_i in seq_along(target_indices)) {
      row <- target_indices[[cell_i]]
      sequence_length <- sum(result$inputs$attention_mask[row, ])
      input_ids_matrix <- result$inputs$gene_ids[row, seq_len(sequence_length), drop = FALSE]
      positions <- which(input_ids_matrix[1L, ] == target_id)
      if (!length(positions)) next
      mask <- matrix(1L, nrow = 1L, ncol = sequence_length)
      input_ids <- torch::torch_tensor(input_ids_matrix, dtype = torch::torch_long(),
                                       device = torch::torch_device("cpu"))
      attention <- torch::torch_tensor(mask, dtype = torch::torch_long(),
                                       device = torch::torch_device("cpu"))
      base <- torch::with_no_grad(model(
        input_ids = input_ids, attention_mask = attention,
        return_layers = TRUE, embedding_layer = -1L
      ))
      altered_ids <- input_ids_matrix
      altered_ids[1L, positions] <- mask_id
      altered <- torch::torch_tensor(altered_ids, dtype = torch::torch_long(),
                                     device = torch::torch_device("cpu"))
      changed <- torch::with_no_grad(model(
        input_ids = altered, attention_mask = attention,
        return_layers = TRUE, embedding_layer = -1L
      ))
      base_array <- torch::as_array(base$hidden_states[[hidden_index]]$cpu())
      changed_array <- torch::as_array(changed$hidden_states[[hidden_index]]$cpu())
      base_hidden <- matrix(base_array[1L, , ], nrow = sequence_length,
                            ncol = model$hidden_size)
      changed_hidden <- matrix(changed_array[1L, , ], nrow = sequence_length,
                               ncol = model$hidden_size)
      token_labels <- fmverse_token_names(result, input_ids_matrix[1L, ])
      eligible <- !(seq_len(sequence_length) %in% positions) &
        !token_labels %in% c("<pad>", "<mask>", "<cls>", "<sep>", "<eos>")
      if (any(eligible)) {
        delta <- sqrt(rowSums((changed_hidden - base_hidden)^2))
        gene_effects[[length(gene_effects) + 1L]] <- data.frame(
          gene = token_labels[eligible], predicted_effect = delta[eligible],
          cell = result$inputs$cell_names[[row]], stringsAsFactors = FALSE
        )
      }
      rm(base, changed, input_ids, altered, attention, base_hidden, changed_hidden,
         base_array, changed_array)
    }
    gene_effects <- do.call(rbind, gene_effects)
    if (is.null(gene_effects) || !nrow(gene_effects)) next
    summaries <- aggregate(predicted_effect ~ gene, gene_effects, mean)
    summaries <- summaries[order(summaries$predicted_effect, decreasing = TRUE), , drop = FALSE]
    summaries <- utils::head(summaries, as.integer(top_n))
    summaries$perturbation <- gene
    outputs[[gene]] <- summaries[, c("perturbation", "gene", "predicted_effect")]
  }
  if (!length(outputs)) stop("No requested Geneformer gene tokens were present in selected cells.")
  effects <- do.call(rbind, outputs)
  effects$gene <- factor(effects$gene, levels = unique(effects$gene))
  heatmap <- ggplot2::ggplot(effects, ggplot2::aes(x = perturbation, y = gene,
                                                   fill = predicted_effect)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::scale_fill_viridis_c(name = "Mean hidden-state shift") +
    ggplot2::labs(
      x = "Masked target gene", y = "Context gene",
      title = "Geneformer in-silico token masking",
      subtitle = "Mean contextual embedding distance after replacing target token with <mask>; exploratory"
    ) + ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  fmverse_save_analysis(result, "perturbation", list(
    effects = effects, cells = selected, cells_by_perturbation = cells_by_perturbation,
    interpretation = paste(
      "Mean per-token hidden-state distance after replacing the target token with <mask>.",
      "This is a context perturbation diagnostic, not a causal expression change."
    )
  ), list(perturbation_heatmap = heatmap), list(
    genes = genes, cells = selected, max_cells = max_cells, top_n = top_n,
    seed = seed, method = "replace token with <mask>; mean hidden-state L2 shift"
  ))
}

#' Infer an attention-derived candidate gene association network
#'
#' @param result `fmverse_result` with saved token inputs.
#' @param cells Optional cells used to average attention.
#' @param max_cells Maximum cells to process.
#' @param max_length Attention sequence cap (default 128 to bound memory).
#' @param top_n Number of edges to return.
#' @param seed Random seed for cell subsampling.
#' @return Updated `fmverse_result` with candidate edge table and network ggplot.
#' @export
.fm_grn_legacy <- function(result, cells = NULL, max_cells = 32L, max_length = 128L,
                   top_n = 100L, seed = 1L) {
  result <- fmverse_check_result(result)
  if (!result$model_name %in% c("scGPT", "Geneformer")) {
    stop("`fm_grn()` currently supports attention extraction from scGPT and ",
      "Geneformer only; ", result$model_name, " has no adapted attention interface.",
      call. = FALSE)
  }
  if (is.null(result$inputs)) stop("Token inputs are absent; re-run fmverse() with the current version.")
  selected <- if (is.null(cells)) result$inputs$cell_names else as.character(cells)
  if (any(!selected %in% result$inputs$cell_names)) stop("Some selected cells are absent from saved tokens.")
  if (length(selected) > max_cells) {
    set.seed(as.integer(seed)); selected <- sample(selected, as.integer(max_cells))
  }
  indices <- match(selected, result$inputs$cell_names)
  inputs <- result$inputs
  active_mask <- if (identical(result$model_name, "Geneformer")) {
    inputs$attention_mask[indices, , drop = FALSE]
  } else {
    !inputs$pad_mask[indices, , drop = FALSE]
  }
  forward <- fmverse_forward_inputs(result, indices, return_attention = TRUE,
                                    max_length = max_length)
  active_length <- forward$sequence_length
  trim <- seq_len(active_length)
  model <- forward$model
  result$model <- model
  attention <- torch::as_array(forward$output$attention$cpu())
  ids <- inputs$gene_ids[indices, trim, drop = FALSE]
  token_names <- matrix(fmverse_token_names(result, as.integer(ids)),
                        nrow = nrow(ids), ncol = ncol(ids))
  edge_sum <- new.env(parent = emptyenv())
  edge_n <- new.env(parent = emptyenv())
  for (cell_index in seq_len(nrow(ids))) {
    valid <- which(active_mask[cell_index, trim] &
                     !token_names[cell_index, ] %in% c("<pad>", "<cls>", "<mask>", "<sep>", "<eos>"))
    if (length(valid) < 2L) next
    local_attention <- attention[cell_index, , , drop = TRUE]
    for (i in valid) for (j in valid) {
      source <- token_names[cell_index, i]
      target <- token_names[cell_index, j]
      if (is.na(source) || is.na(target) || source == target ||
          source %in% c("<pad>", "<cls>") || target %in% c("<pad>", "<cls>")) next
      key <- paste(source, target, sep = "\r")
      previous <- if (exists(key, edge_sum, inherits = FALSE)) get(key, edge_sum) else 0
      previous_n <- if (exists(key, edge_n, inherits = FALSE)) get(key, edge_n) else 0
      assign(key, previous + local_attention[i, j], edge_sum)
      assign(key, previous_n + 1, edge_n)
    }
  }
  keys <- ls(edge_sum, all.names = TRUE)
  if (!length(keys)) stop("No attention edges were available; increase `max_length` or select more cells.")
  edges <- data.frame(
    source = sub("\r.*$", "", keys), target = sub("^.*\r", "", keys),
    weight = vapply(keys, function(key) get(key, edge_sum) / get(key, edge_n), numeric(1)),
    stringsAsFactors = FALSE
  )
  edges <- edges[order(edges$weight, decreasing = TRUE), , drop = FALSE]
  edges <- utils::head(edges, as.integer(top_n))
  # Preserve the full ranked edge table, while limiting the drawn network to
  # its strongest edges and labeling only the more connected genes.
  plot_edges <- utils::head(edges, min(30L, nrow(edges)))
  nodes <- unique(c(plot_edges$source, plot_edges$target))
  angles <- seq(0, 2 * pi, length.out = length(nodes) + 1L)[seq_along(nodes)]
  degree <- table(c(plot_edges$source, plot_edges$target))
  node_data <- data.frame(gene = nodes, x = cos(angles), y = sin(angles),
                          degree = as.integer(degree[nodes]))
  label_nodes <- utils::head(node_data$gene[order(node_data$degree, decreasing = TRUE)],
                             min(18L, nrow(node_data)))
  node_data$label <- ifelse(node_data$gene %in% label_nodes, node_data$gene, "")
  plot_edges$x <- node_data$x[match(plot_edges$source, node_data$gene)]
  plot_edges$y <- node_data$y[match(plot_edges$source, node_data$gene)]
  plot_edges$xend <- node_data$x[match(plot_edges$target, node_data$gene)]
  plot_edges$yend <- node_data$y[match(plot_edges$target, node_data$gene)]
  plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = plot_edges, ggplot2::aes(x, y, xend = xend, yend = yend,
      linewidth = weight), alpha = 0.35, colour = "steelblue") +
    ggplot2::geom_point(data = node_data, ggplot2::aes(x, y, size = degree), colour = "firebrick") +
    ggplot2::geom_text(data = node_data[node_data$label != "", , drop = FALSE],
                       ggplot2::aes(x, y, label = label),
                       nudge_y = 0.04, size = 2.4, check_overlap = TRUE) +
    ggplot2::coord_equal() + ggplot2::scale_linewidth(range = c(0.2, 1.2)) +
    ggplot2::scale_size(range = c(1.4, 3.4), guide = "none") +
    ggplot2::labs(title = paste(result$model_name, "attention-derived candidate associations"),
                  subtitle = paste0("Top ", nrow(plot_edges), " of ", nrow(edges),
                                    " edges; exploratory, not causal regulation"),
                  linewidth = "Attention") +
    ggplot2::theme_void()
  fmverse_save_analysis(result, "grn", list(
    edges = edges[, c("source", "target", "weight")], cells = selected,
    interpretation = "Attention-derived candidate associations are not causal regulatory interactions."
  ), list(grn_network = plot), list(cells = selected, max_cells = max_cells,
                                    max_length = active_length, top_n = top_n, seed = seed))
}

fmverse_binary_auc <- function(truth, score) {
  positive <- truth == 1L
  n_pos <- sum(positive); n_neg <- sum(!positive)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  (sum(rank(score, ties.method = "average")[positive]) - n_pos * (n_pos + 1) / 2) /
    (n_pos * n_neg)
}

#' Predict a drug-versus-control response from cell embeddings
#'
#' @param result `fmverse_result` with treatment labels in cell metadata.
#' @param treatment Metadata column naming treatment or control.
#' @param drug Drug label to distinguish from controls; NULL uses all non-control labels.
#' @param control_label Value used for untreated/vehicle cells.
#' @param dose Optional numeric dose metadata column for a response plot.
#' @param test_fraction Fraction of cells reserved for a reproducible held-out test.
#' @param seed Train/test split seed.
#' @return Updated `fmverse_result` with test predictions, metrics, ROC and optional dose plot.
#' @export
.fm_drug_response_legacy <- function(result, treatment, drug = NULL, control_label = "control",
                             dose = NULL, test_fraction = 0.25, seed = 1L) {
  result <- fmverse_check_result(result)
  meta <- result$seurat[[]]
  if (!treatment %in% colnames(meta)) stop("Treatment metadata column not found: ", treatment)
  if (!is.null(dose) && !dose %in% colnames(meta)) stop("Dose metadata column not found: ", dose)
  labels <- as.character(meta[[treatment]])
  if (is.null(drug)) {
    keep <- !is.na(labels)
    response <- labels != control_label
  } else {
    keep <- !is.na(labels) & labels %in% c(control_label, drug)
    response <- labels == drug
  }
  if (sum(keep & response) < 5L || sum(keep & !response) < 5L) {
    stop("Need at least five treated and five control cells after filtering.")
  }
  cells <- rownames(meta)[keep]
  y <- as.integer(response[keep])
  x <- fmverse_get_embedding(result, "scGPT")[cells, , drop = FALSE]
  set.seed(as.integer(seed))
  split <- unlist(lapply(0:1, function(class) {
    ix <- which(y == class)
    sample(ix, max(1L, min(length(ix) - 1L, floor(length(ix) * test_fraction))))
  }))
  train <- setdiff(seq_along(y), split)
  if (length(unique(y[train])) < 2L || length(unique(y[split])) < 2L) {
    stop("Train and test partitions must each contain treated and control cells.")
  }
  ncomp <- max(1L, min(20L, ncol(x), length(train) - 1L))
  pca <- stats::prcomp(x[train, , drop = FALSE], center = TRUE, scale. = FALSE,
                       rank. = ncomp)
  train_pc <- pca$x[, seq_len(ncomp), drop = FALSE]
  test_pc <- predict(pca, newdata = x[split, , drop = FALSE])[, seq_len(ncomp), drop = FALSE]
  colnames(train_pc) <- colnames(test_pc) <- paste0("PC", seq_len(ncol(train_pc)))
  train_data <- data.frame(response = y[train], train_pc, check.names = FALSE)
  names(train_data)[-1L] <- paste0("PC", seq_len(ncomp))
  fit <- stats::glm(response ~ ., data = train_data, family = stats::binomial())
  probability <- as.numeric(stats::predict(fit, newdata = data.frame(test_pc,
      check.names = FALSE), type = "response"))
  truth <- y[split]
  auc <- fmverse_binary_auc(truth, probability)
  order_ix <- order(probability, decreasing = TRUE)
  precision <- cumsum(truth[order_ix] == 1L) / seq_along(order_ix)
  auprc <- mean(precision[truth[order_ix] == 1L])
  roc <- data.frame(threshold = c(Inf, sort(unique(probability), decreasing = TRUE), -Inf))
  roc$sensitivity <- vapply(roc$threshold, function(t) mean(probability[truth == 1L] >= t), numeric(1))
  roc$specificity <- vapply(roc$threshold, function(t) mean(probability[truth == 0L] < t), numeric(1))
  roc_plot <- ggplot2::ggplot(roc, ggplot2::aes(x = 1 - specificity, y = sensitivity)) +
    ggplot2::geom_path(linewidth = 0.8, colour = "steelblue") +
    ggplot2::geom_abline(linetype = 2, colour = "grey60") +
    ggplot2::coord_equal() + ggplot2::labs(title = "Held-out drug response classification",
      subtitle = sprintf("AUC = %.3f; AUPRC = %.3f", auc, auprc), x = "False positive rate",
      y = "True positive rate") + ggplot2::theme_classic()
  predictions <- data.frame(cell = cells[split], truth = truth,
                            probability = probability, stringsAsFactors = FALSE)
  plots <- list(drug_response_roc = roc_plot)
  if (!is.null(dose)) {
    dose_data <- data.frame(dose = suppressWarnings(as.numeric(meta[cells[split], dose, drop = TRUE])),
                            predicted_response = probability,
                            treatment = labels[match(cells[split], rownames(meta))])
    dose_data <- dose_data[is.finite(dose_data$dose), , drop = FALSE]
    if (nrow(dose_data)) {
      dose_span <- diff(range(dose_data$dose))
      error_width <- if (is.finite(dose_span) && dose_span > 0) dose_span * 0.03 else 0.1
      dose_plot <- ggplot2::ggplot(dose_data, ggplot2::aes(dose, predicted_response,
                                                          colour = treatment)) +
        ggplot2::stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
        ggplot2::stat_summary(fun = mean, geom = "point", size = 2) +
        ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar",
                              width = error_width) +
        ggplot2::labs(title = "Predicted response across observed doses",
                      subtitle = "Held-out dose means; SE shown when multiple cells are available",
                      x = "Dose", y = "Predicted treatment probability",
                      colour = "Treatment") + ggplot2::theme_classic()
      plots$drug_response_dose <- dose_plot
    }
  }
  fmverse_save_analysis(result, "drug_response", list(
    predictions = predictions, metrics = c(auc = auc, auprc = auprc),
    drug = drug, control_label = control_label,
    interpretation = "Performance is estimated for held-out cells from observed treatment conditions; unseen-drug zero-shot prediction is not tested."
  ), plots, list(treatment = treatment, drug = drug, control_label = control_label,
                 dose = dose, test_fraction = test_fraction, seed = seed, n_pcs = ncomp))
}
