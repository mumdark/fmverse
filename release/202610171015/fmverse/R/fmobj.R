# Persistent S4 container. Model computations remain in the legacy kernels.
methods::setClass("fmobj", slots = c(seurat = "Seurat", runs = "list",
  active_run = "character", history = "list", version = "character"))

methods::setValidity("fmobj", function(object) {
  if (length(object@version) != 1L || object@version != "1") return("Unsupported fmobj schema version.")
  if (length(object@active_run) != 1L) return("active_run must be one string.")
  ids <- names(object@runs)
  if (length(object@runs) && (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))) {
    return("Runs must have unique nonempty IDs.")
  }
  if (length(ids) && !object@active_run %in% ids) return("active_run does not identify a run.")
  if (!length(ids) && nzchar(object@active_run)) return("An empty object cannot have an active run.")
  for (type in c("metadata", "reductions", "graphs")) {
    owned <- unlist(lapply(object@runs, function(x) unname(x$artifacts[[type]])), use.names = FALSE)
    if (anyDuplicated(owned)) return(paste("Runs share an owned", type, "artifact."))
  }
  cells <- colnames(object@seurat)
  if (!all(vapply(SeuratObject::Assays(object@seurat), function(name) {
    inherits(object@seurat[[name]], "Assay5")
  }, logical(1)))) return("Shared Seurat assays must use the v5 representation.")
  for (id in ids) {
    run <- object@runs[[id]]
    if (!all(c("model", "metadata", "cells", "reduction", "artifacts", "analyses") %in% names(run))) {
      return(paste("Incomplete run:", id))
    }
    if (anyDuplicated(run$cells) || !all(run$cells %in% cells)) return(paste("Invalid cell IDs:", id))
    if (!run$reduction %in% SeuratObject::Reductions(object@seurat)) return(paste("Missing embedding:", id))
    x <- SeuratObject::Embeddings(object@seurat[[run$reduction]])
    if (!identical(rownames(x), run$cells)) return(paste("Embedding row order differs from run cells:", id))
    for (type in c("metadata", "reductions", "graphs")) {
      available <- switch(type, metadata = colnames(object@seurat[[]]),
        reductions = SeuratObject::Reductions(object@seurat), graphs = names(object@seurat@graphs))
      if (!all(unname(run$artifacts[[type]]) %in% available)) return(paste("Missing", type, "artifact:", id))
    }
  }
  TRUE
})

#' Construct a shared-data foundation-model object
#' @export
fmobj <- function(object) {
  if (methods::is(object, "fmobj")) return(object)
  if (inherits(object, "fmverse_result")) return(as_fmobj(object))
  methods::new("fmobj", seurat = fmverse_to_seurat_v5(object), runs = list(),
    active_run = "", history = list(), version = "1")
}

methods::setMethod("show", "fmobj", function(object) {
  cat(sprintf("fmobj: %d cells, %d features, %d model runs\n",
    ncol(object@seurat), nrow(object@seurat), length(object@runs)))
  if (length(object@runs)) {
    cat("Active run:", object@active_run, "\n")
    print(fm_runs(object), row.names = FALSE)
  }
})

.fm_id <- function(object, run_id = NULL) {
  if (!methods::is(object, "fmobj")) stop("Expected an fmobj object.", call. = FALSE)
  if (is.null(run_id)) run_id <- object@active_run
  if (length(run_id) != 1L || is.na(run_id) || !nzchar(run_id) || !run_id %in% names(object@runs)) {
    stop("Unknown run_id; inspect fm_runs(object) and select a run.", call. = FALSE)
  }
  run_id
}

#' List model runs
#' @export
fm_runs <- function(object) {
  object <- as_fmobj(object)
  if (!length(object@runs)) return(data.frame(run_id = character(), model = character(),
    checkpoint = character(), cells = integer(), dimensions = integer(), active = logical()))
  do.call(rbind, lapply(names(object@runs), function(id) {
    run <- object@runs[[id]]
    data.frame(run_id = id, model = run$model,
      checkpoint = if (is.null(run$metadata$checkpoint_name)) NA_character_ else run$metadata$checkpoint_name,
      cells = length(run$cells), dimensions = ncol(SeuratObject::Embeddings(object@seurat[[run$reduction]])),
      active = identical(id, object@active_run), row.names = NULL)
  }))
}

#' Extract the shared Seurat object
#' @export
fm_seurat <- function(object) as_fmobj(object)@seurat

#' Extract a run record
#' @export
fm_run <- function(object, run_id = NULL) {
  object <- as_fmobj(object)
  object@runs[[.fm_id(object, run_id)]]
}

#' Set the active run without changing its results
#' @export
fm_use <- function(object, run_id) {
  object <- as_fmobj(object)
  object@active_run <- .fm_id(object, run_id)
  object@history[[length(object@history) + 1L]] <- list(task = "use", run_id = run_id,
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE))
  object
}

.fm_fresh <- function(record, task, allow_stale = FALSE) {
  if (is.null(record)) stop("Analysis not found: ", task, ". Run that analysis first.", call. = FALSE)
  if (isTRUE(record$stale) && !allow_stale) stop("Analysis `", task,
    "` is stale because an upstream result changed; rerun it.", call. = FALSE)
  record
}

#' Extract a current analysis, including its parameters and editable plots
#' @export
fm_analysis <- function(object, run_id = NULL, task, allow_stale = FALSE) {
  run <- fm_run(object, run_id)
  .fm_fresh(run$analyses[[task]], task, allow_stale)
}

#' Extract an editable ggplot from a run analysis
#' @export
fm_plot <- function(object, run_id = NULL, task, plot = NULL, allow_stale = FALSE) {
  record <- fm_analysis(object, run_id, task, allow_stale)
  plots <- record$plots
  if (is.null(plot) && length(plots) == 1L) return(plots[[1L]])
  if (!is.null(plot) && plot %in% names(plots)) return(plots[[plot]])
  prefix <- if (task == "reconstruct") "reconstruction" else task
  key <- paste0(prefix, "_", plot)
  if (!is.null(plot) && key %in% names(plots)) return(plots[[key]])
  stop("Select a plot from: ", paste(names(plots), collapse = ", "), call. = FALSE)
}

.fm_embedding_name <- function(object, id, embedding = "original", reduction = NULL) {
  run <- object@runs[[id]]
  if (length(embedding) != 1L || !embedding %in% c("original", "integrated")) {
    stop("embedding must be 'original' or 'integrated'.", call. = FALSE)
  }
  if (!is.null(reduction)) {
    original_aliases <- c(run$reduction, run$model)
    integrated <- run$analyses$integration$reduction
    integrated_alias <- names(run$artifacts$reductions)[match(integrated, unname(run$artifacts$reductions))]
    if (reduction %in% original_aliases) return(run$reduction)
    if (reduction %in% c("fm_harmony", integrated, integrated_alias)) embedding <- "integrated" else {
      stop("reduction does not belong to the selected run; choose the correct run_id.", call. = FALSE)
    }
  }
  if (embedding == "original") return(run$reduction)
  record <- .fm_fresh(run$analyses$integration, "integration")
  record$reduction
}

#' Extract an original or integrated embedding in the run's cell order
#' @export
fm_embeddings <- function(object, run_id = NULL, embedding = "original") {
  object <- as_fmobj(object)
  id <- .fm_id(object, run_id)
  name <- .fm_embedding_name(object, id, embedding)
  x <- SeuratObject::Embeddings(object@seurat[[name]])[object@runs[[id]]$cells, , drop = FALSE]
  if (embedding == "original") colnames(x) <- object@runs[[id]]$embedding_columns
  x
}

.fm_translate <- function(x, mapping) {
  if (is.character(x) && length(mapping)) {
    hit <- match(x, names(mapping))
    x[!is.na(hit)] <- unname(mapping[hit[!is.na(hit)]])
  } else if (is.list(x) && !is.object(x)) {
    x <- lapply(x, .fm_translate, mapping = mapping)
  }
  x
}

.fm_reverse <- function(mapping) {
  if (length(mapping)) stats::setNames(names(mapping), unname(mapping)) else character()
}

.fm_run_id <- function(object, result, run_id) {
  if (!is.null(run_id)) {
    if (!is.character(run_id) || length(run_id) != 1L || is.na(run_id) || !nzchar(run_id)) {
      stop("run_id must be one nonempty string.", call. = FALSE)
    }
    return(run_id)
  }
  checkpoint <- result$metadata$checkpoint_name
  if (is.null(checkpoint) || length(checkpoint) != 1L || is.na(checkpoint) || !nzchar(checkpoint)) checkpoint <- "local"
  base <- paste(result$model_name, checkpoint, sep = "_")
  id <- base; i <- 1L
  while (id %in% names(object@runs)) { i <- i + 1L; id <- paste0(base, "_", i) }
  id
}

.fm_drop_artifacts <- function(seurat, artifacts) {
  for (type in names(artifacts)) for (name in unname(artifacts[[type]])) {
    available <- switch(type, metadata = colnames(seurat[[]]),
      reductions = SeuratObject::Reductions(seurat), graphs = names(seurat@graphs))
    if (name %in% available) seurat[[name]] <- NULL
  }
  seurat
}

# Build a short-lived legacy view; no algorithm has to know about the registry.
.fm_view <- function(object, id) {
  run <- object@runs[[id]]
  seurat <- object@seurat
  for (other in object@runs) seurat <- .fm_drop_artifacts(seurat, other$artifacts)
  if (!identical(colnames(seurat), run$cells)) seurat <- subset(seurat, cells = run$cells)
  for (local in names(run$artifacts$metadata)) {
    physical <- run$artifacts$metadata[[local]]
    value <- object@seurat[[]][run$cells, physical, drop = FALSE]
    names(value) <- local
    seurat <- SeuratObject::AddMetaData(seurat, value)
  }
  for (local in names(run$artifacts$reductions)) {
    dr <- object@seurat[[run$artifacts$reductions[[local]]]]
    SeuratObject::Key(dr) <- paste0(gsub("[^A-Za-z0-9]", "", local), "_")
    seurat[[local]] <- dr
  }
  for (local in names(run$artifacts$graphs)) seurat[[local]] <- object@seurat[[run$artifacts$graphs[[local]]]]
  mapping <- do.call(c, unname(lapply(run$artifacts, .fm_reverse)))
  analyses <- lapply(run$analyses, .fm_translate, mapping = mapping)
  analyses <- analyses[!vapply(analyses, function(x) isTRUE(x$stale), logical(1))]
  plots <- unlist(lapply(analyses, function(x) x$plots), recursive = FALSE)
  if (length(plots)) names(plots) <- unlist(lapply(analyses, function(x) names(x$plots)), use.names = FALSE)
  structure(list(model_name = run$model, embeddings = fm_embeddings(object, id),
    seurat = seurat, model = NULL, audit = run$audit, inputs = run$inputs,
    metadata = run$metadata, analyses = analyses, plots = plots, history = list()),
    class = c("fmverse_result", "list"))
}

.fm_artifact_name <- function(prefix, type, local) {
  paste0(prefix, "_", substr(type, 1L, 1L), "_", digest::digest(local, algo = "crc32", serialize = FALSE))
}

.fm_legacy_artifacts <- function(result, type) {
  a <- result$analyses
  switch(type,
    metadata = unique(c(a$cluster$cluster_column,
      if (!is.null(a$cluster)) c("seurat_clusters", grep("^fm_tmp_", colnames(result$seurat[[]]), value = TRUE)),
      a$annotation$prediction_column, a$annotation$confidence_column, colnames(a$trajectory$pseudotime))),
    reductions = unique(c(a$integration$reduction,
      if (!is.null(a$cluster)) "umap.fm_cluster", if (!is.null(a$annotation)) "umap.fm_annotation",
      if (!is.null(a$trajectory)) "umap.fm_trajectory",
      if (!is.null(a$integration)) c("umap.fm_before_integration", "umap.fm_after_integration"))),
    graphs = if (!is.null(a$cluster)) grep("^fm_(leiden|louvain)_(nn|snn)$", names(result$seurat@graphs), value = TRUE) else character())
}

# Merge only run-owned artifacts; other runs and the shared raw data survive.
.fm_store_view <- function(object, id, result, before, migration = FALSE) {
  run <- object@runs[[id]]
  source <- result$seurat
  for (type in c("metadata", "reductions", "graphs")) {
    names_of <- function(s) switch(type, metadata = colnames(s[[]]),
      reductions = SeuratObject::Reductions(s), graphs = names(s@graphs))
    src <- names_of(source); prev <- names_of(before)
    missing <- setdiff(names(run$artifacts[[type]]), src)
    if (length(missing)) {
      for (physical in unname(run$artifacts[[type]][missing])) object@seurat[[physical]] <- NULL
      run$artifacts[[type]] <- run$artifacts[[type]][!names(run$artifacts[[type]]) %in% missing]
    }
    locals <- union(names(run$artifacts[[type]]), setdiff(src, prev))
    common <- intersect(src, prev)
    changed <- common[!vapply(common, function(n) {
      if (type == "metadata") identical(source[[n, drop = TRUE]], before[[n, drop = TRUE]])
      else identical(source[[n]], before[[n]])
    }, logical(1))]
    locals <- union(locals, changed)
    if (type == "reductions") locals <- union(locals, run$model)
    if (migration) locals <- union(locals, .fm_legacy_artifacts(result, type))
    locals <- intersect(locals, src)
    for (local in locals) {
      physical <- unname(run$artifacts[[type]][local])
      if (!length(physical) || is.na(physical)) physical <- .fm_artifact_name(run$prefix, type, local)
      if (type == "reductions" && identical(local, run$model)) physical <- run$reduction
      run$artifacts[[type]][local] <- physical
      if (type == "metadata") {
        value <- source[[]][, local, drop = FALSE]; names(value) <- physical
        object@seurat <- SeuratObject::AddMetaData(object@seurat, value)
      } else {
        value <- source[[local]]
        if (type == "reductions") SeuratObject::Key(value) <- paste0(toupper(gsub("[^A-Za-z0-9]", "", physical)), "_")
        object@seurat[[physical]] <- value
      }
    }
  }
  if (identical(colnames(source), colnames(object@seurat))) {
    # Preserve normalization produced by existing Seurat algorithms, once.
    for (assay in SeuratObject::Assays(source)) object@seurat[[assay]] <- source[[assay]]
  }
  object@runs[[id]] <- run
  object
}

.fm_add_legacy <- function(object, result, run_id = NULL, overwrite = FALSE, migration = FALSE) {
  id <- .fm_run_id(object, result, run_id)
  if (id %in% names(object@runs)) {
    if (!isTRUE(overwrite)) stop("run_id already exists; choose another ID or set overwrite = TRUE.", call. = FALSE)
    prefix <- object@runs[[id]]$prefix
    object@seurat <- .fm_drop_artifacts(object@seurat, object@runs[[id]]$artifacts)
  } else {
    occupied <- vapply(object@runs, function(x) x$prefix, character(1))
    number <- 1L
    repeat {
      prefix <- sprintf("fmR%04d", number)
      used_names <- c(colnames(object@seurat[[]]), SeuratObject::Reductions(object@seurat), names(object@seurat@graphs))
      if (!prefix %in% occupied && !prefix %in% used_names && !any(startsWith(used_names, paste0(prefix, "_")))) break
      number <- number + 1L
    }
  }
  run <- list(model = result$model_name, metadata = result$metadata,
    parameters = result$metadata, cells = rownames(result$embeddings),
    reduction = prefix, prefix = prefix, embedding_columns = colnames(result$embeddings),
    audit = result$audit, preprocessing = result$metadata$preprocessing,
    inputs = result$inputs, artifacts = list(metadata = character(), reductions = character(), graphs = character()),
    analyses = list())
  run$artifacts$reductions[result$model_name] <- prefix
  object@runs[[id]] <- run
  before <- object@seurat
  object <- .fm_store_view(object, id, result, before, migration)
  if (migration) {
    locals <- object@runs[[id]]$artifacts
    protected <- unlist(lapply(object@runs, function(x) unlist(x$artifacts, use.names = FALSE)), use.names = FALSE)
    for (type in names(locals)) for (name in setdiff(names(locals[[type]]), protected)) {
      if (name %in% c(colnames(object@seurat[[]]), SeuratObject::Reductions(object@seurat), names(object@seurat@graphs))) {
        object@seurat[[name]] <- NULL
      }
    }
  }
  mapping <- do.call(c, unname(object@runs[[id]]$artifacts))
  for (task in names(result$analyses)) {
    record <- .fm_translate(result$analyses[[task]], mapping)
    record$embeddings <- NULL
    record$stale <- FALSE; record$revision <- 1L; record$depends_on <- character()
    if (task == "enrichment" && !is.null(result$analyses$markers)) record$depends_on <- "markers"
    if (task %in% c("markers", "trajectory") && !is.null(result$analyses$cluster) &&
        identical(result$analyses[[task]]$group_by, result$analyses$cluster$cluster_column)) record$depends_on <- "cluster"
    prefix_plot <- if (task == "reconstruct") "reconstruction_" else paste0(task, "_")
    record$plots <- lapply(result$plots[startsWith(names(result$plots), prefix_plot)], fmverse_detach_plot)
    record$parameters <- list()
    object@runs[[id]]$analyses[[task]] <- record
  }
  object@active_run <- id
  if (migration && length(result$history)) {
    object@history <- c(object@history, lapply(result$history, function(entry) { entry$run_id <- id; entry }))
  }
  object@history[[length(object@history) + 1L]] <- list(task = if (migration) "migrate" else "inference",
    run_id = id, timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE), parameters = result$metadata)
  methods::validObject(object)
  object
}

#' Convert an old single-run result without recomputing inference
#' @export
as_fmobj <- function(object, run_id = NULL) {
  if (methods::is(object, "fmobj")) return(object)
  if (inherits(object, "Seurat")) return(fmobj(object))
  if (!inherits(object, "fmverse_result")) stop("Expected fmobj, Seurat or legacy fmverse_result.", call. = FALSE)
  object$seurat <- fmverse_to_seurat_v5(object$seurat)
  if (!object$model_name %in% SeuratObject::Reductions(object$seurat)) {
    object$seurat[[object$model_name]] <- SeuratObject::CreateDimReducObject(
      embeddings = object$embeddings, key = "FMBASELINE_", assay = SeuratObject::DefaultAssay(object$seurat))
  }
  .fm_add_legacy(fmobj(object$seurat), object, run_id, migration = TRUE)
}
