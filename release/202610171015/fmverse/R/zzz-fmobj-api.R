.fm_inference_model_dir_api <- function(object, model_dir, dots, run_id = NULL,
                                        overwrite = FALSE, expected_model = NULL) {
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    stop("All inference options after `model_dir` must be named.", call. = FALSE)
  }
  forbidden <- intersect(names(dots), c("model", "checkpoint", "checkpoint_dir",
    "vocab_file", "gene_index_file", "config_file", "assets_dir"))
  if (length(forbidden)) {
    stop("Model selection and resource paths are read from `fmverse.json`; remove `",
         paste(forbidden, collapse = "`, `"), "`.", call. = FALSE)
  }
  duplicated <- unique(names(dots)[duplicated(names(dots))])
  if (length(duplicated)) {
    stop("Inference options may only be supplied once: ",
         paste(duplicated, collapse = ", "), call. = FALSE)
  }
  allowed <- setdiff(names(formals(.fmverse_infer_legacy)),
                     c("object", "model", "checkpoint", "checkpoint_dir"))
  unknown <- setdiff(names(dots), allowed)
  if (length(unknown)) {
    stop("Unknown inference option(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  info <- .fm_model_manifest(model_dir)
  if (!is.null(expected_model) && !identical(info$model, expected_model)) {
    stop("This shortcut requires a ", expected_model,
         " model directory; selected directory is marked `", info$model, "`.",
         call. = FALSE)
  }
  object <- as_fmobj(object)
  if (!is.null(run_id) && (!is.character(run_id) || length(run_id) != 1L ||
      is.na(run_id) || !nzchar(run_id))) {
    stop("run_id must be one nonempty string.", call. = FALSE)
  }
  if (!is.null(run_id) && run_id %in% names(object@runs) && !overwrite) {
    stop("run_id already exists; choose another ID or set overwrite = TRUE.",
         call. = FALSE)
  }
  source <- object@seurat
  if (!is.null(run_id) && run_id %in% names(object@runs) && overwrite) {
    source <- .fm_drop_artifacts(source, object@runs[[run_id]]$artifacts)
  }

  if (is.null(dots$batch_size)) {
    dots$batch_size <- switch(info$model, scGPT = 16L, Geneformer = 8L, 1L)
  }
  if (identical(info$model, "scGPT") && is.null(dots$vocab_file)) {
    dots$vocab_file <- info$files$vocab
    dots$config_file <- info$files$config
  }
  if (identical(info$model, "scFoundation")) {
    dots$config_file <- info$files$config
  }
  if (info$model %in% c("scBERT", "scFoundation")) {
    dots$gene_index_file <- info$files$gene_index
  }
  if (identical(info$model, "UCE") && is.null(dots$species)) {
    dots$species <- info$defaults$species
  }
  if (identical(info$model, "UCE")) dots$assets_dir <- info$files$assets_dir
  checkpoint_path <- switch(info$model,
    scGPT = info$files$weights,
    Geneformer = info$files$weights,
    scBERT = info$files$weights,
    scFoundation = info$files$weights,
    UCE = info$files$weights
  )
  if (identical(info$model, "Geneformer")) dots$checkpoint_dir <- info$root
  kernel_args <- c(list(object = source, model = info$model,
                        checkpoint = checkpoint_path), dots)
  result <- do.call(.fmverse_infer_legacy, kernel_args)
  result$metadata$checkpoint_name <- info$checkpoint
  result$metadata$model_dir <- info$root
  result$metadata$manifest <- file.path(info$root, "fmverse.json")
  .fm_add_legacy(object, result, run_id, overwrite)
}

#' Run a registered R-native model from a self-contained model directory
#'
#' The selected folder must contain `fmverse.json`, which declares its model,
#' checkpoint, and relative paths to all required files. The Seurat object (or
#' existing `fmobj`) remains the first argument; `model` and `checkpoint` are
#' inferred from the folder.
#'
#' @param object A Seurat v4/v5 object or an existing `fmobj`.
#' @param model_dir Directory containing `fmverse.json` and model resources.
#' @param ... Named model options such as `batch_size`, `cells`, `assay`,
#'   `device`, or model-specific preprocessing options.
#' @param run_id Optional name for this run. Defaults to model plus checkpoint.
#' @param overwrite Replace an existing run with the same explicit ID.
#' @return An `fmobj` containing the new model run.
#' @export
fmverse <- function(object, model_dir, ..., run_id = NULL, overwrite = FALSE) {
  .fm_inference_model_dir_api(object, model_dir, list(...), run_id, overwrite)
}

#' @rdname fmverse
#' @export
fmverse_geneformer <- function(object, model_dir, ..., run_id = NULL,
                               overwrite = FALSE) {
  .fm_inference_model_dir_api(object, model_dir, list(...), run_id, overwrite,
                              expected_model = "Geneformer")
}

#' @rdname fmverse
#' @export
fmverse_scbert <- function(object, model_dir, ..., run_id = NULL,
                           overwrite = FALSE) {
  .fm_inference_model_dir_api(object, model_dir, list(...), run_id, overwrite,
                              expected_model = "scBERT")
}

#' @rdname fmverse
#' @export
fmverse_scfoundation <- function(object, model_dir, ..., run_id = NULL,
                                 overwrite = FALSE) {
  .fm_inference_model_dir_api(object, model_dir, list(...), run_id, overwrite,
                              expected_model = "scFoundation")
}

#' @rdname fmverse
#' @export
fmverse_uce <- function(object, model_dir, ..., run_id = NULL,
                        overwrite = FALSE) {
  .fm_inference_model_dir_api(object, model_dir, list(...), run_id, overwrite,
                              expected_model = "UCE")
}

.fm_invalidate <- function(analyses, task) {
  changed <- task
  repeat {
    dependent <- names(analyses)[vapply(analyses,
      function(x) any(x$depends_on %in% changed), logical(1))]
    added <- setdiff(dependent, changed)
    if (!length(added)) break
    for (name in added) analyses[[name]]$stale <- TRUE
    changed <- union(changed, added)
  }
  analyses
}

.fm_dispatch <- function(task, kernel, args) {
  object <- as_fmobj(args$result)
  id <- .fm_id(object, args$run_id)
  run <- object@runs[[id]]
  view <- .fm_view(object, id)
  embedding <- if (is.null(args$embedding)) "original" else args$embedding
  dependencies <- character()
  embedding_tasks <- c("cluster", "annotation", "integration", "trajectory", "drug_response")
  if (task %in% embedding_tasks) {
    physical <- .fm_embedding_name(object, id, embedding, args[["reduction"]])
    local <- names(run$artifacts$reductions)[match(physical, unname(run$artifacts$reductions))]
    if (!length(local) || is.na(local)) stop("The selected run has no matching embedding.")
    embedding <- if (identical(physical, run$reduction)) "original" else "integrated"
    if (!identical(physical, run$reduction)) dependencies <- c(dependencies, "integration")
    if ("reduction" %in% names(formals(kernel))) args$reduction <- local
  }
  reverse_meta <- .fm_reverse(run$artifacts$metadata)
  for (name in c("group_by", "batch", "truth")) {
    if (!is.null(args[[name]])) args[[name]] <- .fm_translate(args[[name]], reverse_meta)
  }
  if (task %in% c("markers", "trajectory") && is.null(args$group_by)) args$group_by <- "fm_cluster"
  group <- args$group_by
  if (!is.null(group)) {
    cluster_column <- .fm_translate(run$analyses$cluster$cluster_column, reverse_meta)
    annotation_column <- .fm_translate(run$analyses$annotation$prediction_column, reverse_meta)
    if (identical(group, "fm_cluster") || (!is.null(cluster_column) && identical(group, cluster_column))) {
      cluster <- .fm_fresh(run$analyses$cluster, "cluster")
      args$group_by <- .fm_translate(cluster$cluster_column, reverse_meta)
      dependencies <- c(dependencies, "cluster")
    }
    if (!is.null(annotation_column) && identical(group, annotation_column)) {
      .fm_fresh(run$analyses$annotation, "annotation")
      dependencies <- c(dependencies, "annotation")
    }
  }
  if (task == "enrichment" && is.null(args$markers)) {
    .fm_fresh(run$analyses$markers, "markers")
    dependencies <- c(dependencies, "markers")
  }
  if (task == "annotation") {
    reference <- as_fmobj(args$reference)
    reference_id <- .fm_id(reference, args$reference_run_id)
    reference_run <- reference@runs[[reference_id]]
    reference_physical <- .fm_embedding_name(reference, reference_id, embedding)
    reference_local <- names(reference_run$artifacts$reductions)[match(reference_physical,
      unname(reference_run$artifacts$reductions))]
    ref_view <- .fm_view(reference, reference_id)
    # Both temporary views use the query reduction alias; registry IDs stay private.
    if (!identical(local, reference_local)) ref_view$seurat[[local]] <- ref_view$seurat[[reference_local]]
    args$reference <- ref_view
    ref_meta <- .fm_reverse(reference_run$artifacts$metadata)
    if (!is.null(args$label)) args$label <- .fm_translate(args$label, ref_meta)
  }
  if (task == "integration" && !is.null(args$reduction_name) &&
      args$reduction_name %in% c(run$model, run$reduction)) {
    stop("Integration must write a new reduction; the original embedding is protected.", call. = FALSE)
  }
  # Drop this task's former outputs in the temporary view. The merge removes
  # abandoned physical artifacts, including renamed columns and graph methods.
  if (!is.null(run$analyses[[task]])) {
    old <- view
    reverse <- do.call(c, unname(lapply(run$artifacts, .fm_reverse)))
    old$analyses <- stats::setNames(list(.fm_translate(run$analyses[[task]], reverse)), task)
    protected <- c(run$model, args[["reduction"]], args$group_by, args$batch, args$truth)
    for (type in c("metadata", "reductions", "graphs")) {
      obsolete <- intersect(.fm_legacy_artifacts(old, type), names(run$artifacts[[type]]))
      for (name in setdiff(obsolete, protected)) view$seurat[[name]] <- NULL
    }
  }
  # Existing UMAP helpers cache by name. Remove only this task's display before
  # rerunning so changing dimensions or embedding cannot reuse obsolete coordinates.
  display <- switch(task, cluster = "umap.fm_cluster", annotation = "umap.fm_annotation",
    trajectory = "umap.fm_trajectory", integration = c("umap.fm_before_integration", "umap.fm_after_integration"),
    character())
  for (name in intersect(display, SeuratObject::Reductions(view$seurat))) view$seurat[[name]] <- NULL
  drug_original <- NULL
  if (task == "drug_response") {
    if ("scGPT" %in% SeuratObject::Reductions(view$seurat)) drug_original <- view$seurat[["scGPT"]]
    view$seurat[["scGPT"]] <- view$seurat[[local]]
  }
  before <- view$seurat
  args$result <- view
  args$run_id <- args$reference_run_id <- args$embedding <- NULL
  # A NULL compatibility reduction means use the selected run, never guess a model.
  if ("reduction" %in% names(args) && is.null(args[["reduction"]])) args$reduction <- local
  result <- do.call(kernel, args)
  if (task == "drug_response") {
    if (is.null(drug_original)) result$seurat[["scGPT"]] <- NULL else result$seurat[["scGPT"]] <- drug_original
  }
  object <- .fm_store_view(object, id, result, before)
  record <- result$analyses[[task]]
  mapping <- do.call(c, unname(object@runs[[id]]$artifacts))
  record <- .fm_translate(record, mapping)
  record$embeddings <- NULL
  record$parameters <- .fm_translate(tail(result$history, 1L)[[1L]]$parameters, mapping)
  record$parameters$embedding <- embedding
  record$parameters$run_id <- id
  if (task == "annotation") record$parameters$reference_run_id <- reference_id
  prefix <- if (task == "reconstruct") "reconstruction_" else paste0(task, "_")
  record$plots <- lapply(result$plots[startsWith(names(result$plots), prefix)], fmverse_detach_plot)
  record$revision <- if (is.null(run$analyses[[task]]$revision)) 1L else run$analyses[[task]]$revision + 1L
  record$depends_on <- unique(setdiff(dependencies, task))
  record$stale <- FALSE
  analyses <- .fm_invalidate(object@runs[[id]]$analyses, task)
  analyses[[task]] <- record
  object@runs[[id]]$analyses <- analyses
  object@history[[length(object@history) + 1L]] <- list(task = task, run_id = id,
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE), parameters = record$parameters)
  methods::validObject(object)
  object
}

.fm_downstream_wrapper <- function(task, kernel) {
  wrapper <- function() {
    supplied <- names(as.list(match.call())[-1L])
    args <- mget(supplied, envir = environment(), inherits = FALSE)
    .fm_dispatch(task, kernel, args)
  }
  f <- formals(kernel)
  if ("reduction" %in% names(f)) f["reduction"] <- alist(reduction = NULL)
  if (task %in% c("markers", "trajectory")) f["group_by"] <- alist(group_by = NULL)
  f <- c(f, alist(run_id = NULL))
  if (task %in% c("cluster", "annotation", "integration", "trajectory", "drug_response")) {
    f <- c(f, alist(embedding = "original"))
  }
  if (task == "annotation") f <- c(f, alist(reference_run_id = NULL))
  formals(wrapper) <- f
  wrapper
}

fm_cluster <- .fm_downstream_wrapper("cluster", .fm_cluster_legacy)
fm_annotate <- .fm_downstream_wrapper("annotation", .fm_annotate_legacy)
fm_integrate <- .fm_downstream_wrapper("integration", .fm_integrate_legacy)
fm_markers <- .fm_downstream_wrapper("markers", .fm_markers_legacy)
fm_enrich <- .fm_downstream_wrapper("enrichment", .fm_enrich_legacy)
fm_trajectory <- .fm_downstream_wrapper("trajectory", .fm_trajectory_legacy)
fm_reconstruct <- .fm_downstream_wrapper("reconstruct", .fm_reconstruct_legacy)
fm_perturb <- .fm_downstream_wrapper("perturbation", .fm_perturb_legacy)
fm_grn <- .fm_downstream_wrapper("grn", .fm_grn_legacy)
fm_drug_response <- .fm_downstream_wrapper("drug_response", .fm_drug_response_legacy)
