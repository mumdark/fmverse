# scFoundation gene alignment, count preprocessing, and cell embedding inference.

fmverse_scfoundation_gene_index <- function(gene_index_file = NULL) {
  if (is.null(gene_index_file)) {
    gene_index_file <- system.file("extdata", "scfoundation_gene_index_19264.tsv",
      package = "fmverse")
    if (!nzchar(gene_index_file)) {
      candidate <- file.path("inst", "extdata", "scfoundation_gene_index_19264.tsv")
      if (file.exists(candidate)) gene_index_file <- candidate
    }
  }
  if (!length(gene_index_file) || !nzchar(gene_index_file) || !file.exists(gene_index_file)) {
    stop("scFoundation 19,264-gene index is missing from the package or selected path.",
      call. = FALSE)
  }
  tab <- utils::read.delim(gene_index_file, stringsAsFactors = FALSE,
    check.names = FALSE, quote = "", comment.char = "")
  if (!all(c("gene_name", "index") %in% names(tab))) {
    stop("scFoundation gene index must contain `gene_name` and `index` columns.")
  }
  tab <- tab[order(as.integer(tab$index)), , drop = FALSE]
  genes <- trimws(as.character(tab$gene_name))
  if (length(genes) != 19264L || anyNA(genes) || any(!nzchar(genes)) ||
      anyDuplicated(genes) || !identical(as.integer(tab$index), 0:19263)) {
    stop("scFoundation gene index must have 19,264 unique genes with sequential zero-based indices.")
  }
  genes
}

fmverse_scfoundation_tokenize <- function(object, assay = NULL, cells = NULL,
                                           gene_index_file = NULL,
                                           pre_normalized = "F", tgthighres = "a5") {
  pre_normalized <- match.arg(pre_normalized, c("F", "T"))
  gene_order <- fmverse_scfoundation_gene_index(gene_index_file)
  object <- fmverse_to_seurat_v5(object)
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
  if (!assay %in% SeuratObject::Assays(object)) stop("Assay not found: ", assay)
  if (!"counts" %in% SeuratObject::Layers(object[[assay]])) {
    stop("scFoundation input requires a Seurat `counts` layer.")
  }
  counts <- methods::as(SeuratObject::LayerData(object, assay = assay, layer = "counts"),
    "dgCMatrix")
  if (is.null(cells)) cell_names <- colnames(counts) else {
    cell_names <- as.character(cells)
    if (!length(cell_names) || anyDuplicated(cell_names) || any(!cell_names %in% colnames(counts))) {
      stop("`cells` must contain unique barcodes present in the Seurat object.")
    }
    counts <- counts[, cell_names, drop = FALSE]
  }
  feature_metadata <- object[[assay]][[]]
  symbols <- if ("scGPT_gene_name" %in% colnames(feature_metadata)) {
    as.character(feature_metadata$scGPT_gene_name)
  } else rownames(counts)
  gene_index <- match(symbols, gene_order)
  valid <- which(!is.na(gene_index))
  if (!length(valid)) stop("No Seurat features match the scFoundation gene index.")
  collapse <- Matrix::sparseMatrix(i = gene_index[valid], j = seq_along(valid), x = 1,
    dims = c(length(gene_order), length(valid)))
  aligned <- methods::as(collapse %*% counts[valid, , drop = FALSE], "dgCMatrix")
  rownames(aligned) <- gene_order
  colnames(aligned) <- cell_names
  if (grepl("^[aft][0-9]+(?:[.][0-9]+)?$", tgthighres, ignore.case = TRUE) == FALSE) {
    stop("`tgthighres` must be an official resolution string such as `a5`, `f1`, or `t4`.")
  }
  resolution_mode <- tolower(substr(tgthighres, 1L, 1L))
  resolution_value <- as.numeric(sub("^[aft]", "", tolower(tgthighres)))
  n_cells <- ncol(aligned)
  cell_inputs <- vector("list", n_cells)
  token_counts <- integer(n_cells)
  totals <- numeric(n_cells)
  pointers <- aligned@p
  for (j in seq_len(n_cells)) {
    first <- pointers[[j]] + 1L; last <- pointers[[j + 1L]]
    if (first > last) stop("Cell `", cell_names[[j]], "` has no counts in the scFoundation gene set.")
    offsets <- seq.int(first, last)
    gene_pos <- aligned@i[offsets] + 1L
    expression <- as.numeric(aligned@x[offsets])
    total_count <- sum(expression)
    if (!is.finite(total_count) || total_count <= 0) {
      stop("Cell `", cell_names[[j]], "` has nonpositive scFoundation total counts.")
    }
    totals[[j]] <- total_count
    gene_values <- if (pre_normalized == "F") {
      log1p(expression / total_count * 1e4)
    } else expression
    total_first <- switch(resolution_mode,
      a = log10(total_count) + resolution_value,
      f = log10(total_count * resolution_value),
      t = resolution_value)
    all_values <- c(gene_values, total_first, log10(total_count))
    all_positions <- c(gene_pos - 1L, 19264L, 19265L)
    keep <- all_values > 0
    values <- all_values[keep]
    positions <- all_positions[keep]
    # `gatherData` compacts positive entries in the original fixed-gene order;
    # two appended total-count tokens remain the final sequence elements.
    order_idx <- order(positions, method = "radix")
    values <- values[order_idx]
    positions <- positions[order_idx]
    special <- positions >= 19264L
    values <- c(values[!special], values[special])
    positions <- c(positions[!special], positions[special])
    if (length(values) < 2L || !all(positions[(length(positions) - 1L):length(positions)] >= 19264L)) {
      stop("scFoundation total-count tokens are absent; check `tgthighres` and cell totals.")
    }
    cell_inputs[[j]] <- list(values = values, gene_positions = as.integer(positions))
    token_counts[[j]] <- length(values)
  }
  list(cells = cell_inputs, cell_names = cell_names, gene_order = gene_order,
    assay = assay, matched_features = length(valid), total_features = nrow(counts),
    token_counts = token_counts, total_counts = totals,
    preprocessing = list(method = "scFoundation official singlecell preprocessing",
      input = if (pre_normalized == "F") "raw counts" else "pre-normalized expression",
      normalization = if (pre_normalized == "F") "log1p(counts / mapped_total * 1e4)" else "identity",
      resolution = tgthighres, compacting = "retain positive gene values in fixed order; append two total-count tokens",
      gene_index_file = if (is.null(gene_index_file)) "package bundled scFoundation 19,264-gene index" else
        normalizePath(gene_index_file, winslash = "/", mustWork = TRUE)))
}

fmverse_scfoundation_checkpoint_config <- function(resolved, config_file = NULL) {
  if (!is.null(config_file)) {
    if (!file.exists(config_file)) stop("scFoundation config file does not exist: ", config_file,
      call. = FALSE)
    return(jsonlite::fromJSON(config_file, simplifyVector = FALSE))
  }
  candidates <- file.path(resolved$directory, c("scfoundation_cell_config.json", "config.json"))
  found <- candidates[file.exists(candidates)]
  if (!length(found)) stop("scFoundation model config was not exported beside the flat tensor file. ",
    "Run 00scFoundation/prepare_checkpoint.py on the official MMF checkpoint.", call. = FALSE)
  jsonlite::fromJSON(found[[1L]], simplifyVector = FALSE)
}

#' Embed cells with the R-native scFoundation cell checkpoint.
#' @param object Seurat v4/v5 or `fmobj`.
#' @param model Model option, `scFoundation`.
#' @param checkpoint Official checkpoint name `01B-resolution` or a local file.
#' @param checkpoint_dir Optional checkpoint root.
#' @param gene_index_file Optional fixed 19,264-gene index; the official list is bundled.
#' @param config_file Internal explicit path to the scFoundation config from
#'   the selected folder manifest.
#' @param batch_size Number of cells per progress chunk. The official reference computes
#'   one variable-length cell sequence at a time.
#' @param cells Optional cell barcodes.
#' @param assay Counts assay.
#' @param device Torch device.
#' @param verbose Report inference progress.
#' @param pre_normalized `F` for raw counts or `T` for pre-normalized expression.
#' @param tgthighres Official total-count resolution token, for example `a5`.
#' @param pool_type `all` (four-part 2,048D by default) or `max`.
#' @return An `fmverse_result` suitable for conversion to an `fmobj` run.
#' @export
.fmverse_scfoundation_legacy <- function(object, model = "scFoundation",
    checkpoint = "01B-resolution", checkpoint_dir = NULL, gene_index_file = NULL,
    config_file = NULL,
    batch_size = 1L, cells = NULL, assay = NULL, device = "cpu",
    verbose = interactive(), pre_normalized = "F", tgthighres = "a5",
    pool_type = "all") {
  if (length(batch_size) != 1L || is.na(batch_size) || batch_size < 1L ||
      batch_size != as.integer(batch_size)) stop("`batch_size` must be a positive integer.")
  pool_type <- match.arg(pool_type, c("all", "max"))
  resolved <- fmverse_resolve_scfoundation_checkpoint(checkpoint, checkpoint_dir)
  config <- fmverse_scfoundation_checkpoint_config(resolved, config_file)
  tokenized <- fmverse_scfoundation_tokenize(object, assay, cells, gene_index_file,
    pre_normalized = pre_normalized, tgthighres = tgthighres)
  model_module <- fmverse_scfoundation_build_model(config)
  loaded <- fmverse_load_scfoundation_weights(model_module, resolved$path)
  model_module <- loaded$model
  torch_device <- torch::torch_device(device)
  model_module$to(device = torch_device); model_module$eval()
  n_cells <- length(tokenized$cell_names)
  hidden_size <- as.integer(config$encoder$hidden_dim)
  output_dim <- if (pool_type == "all") 4L * hidden_size else hidden_size
  embeddings <- matrix(NA_real_, nrow = n_cells, ncol = output_dim)
  starts <- seq.int(1L, n_cells, by = as.integer(batch_size))
  progress_step <- max(1L, ceiling(length(starts) / 20L))
  for (batch_no in seq_along(starts)) {
    first <- starts[[batch_no]]
    index <- first:min(first + as.integer(batch_size) - 1L, n_cells)
    for (cell_index in index) {
      cell <- tokenized$cells[[cell_index]]
      value_tensor <- torch::torch_tensor(matrix(as.numeric(cell$values), nrow = 1L),
        dtype = torch::torch_float(), device = torch_device)
      position_tensor <- torch::torch_tensor(matrix(as.integer(cell$gene_positions), nrow = 1L),
        dtype = torch::torch_long(), device = torch_device)
      hidden <- torch::with_no_grad(model_module(values = value_tensor,
        gene_positions = position_tensor))
      seq_len <- hidden$size(2L)
      if (seq_len < 3L) stop("scFoundation produced fewer than two total tokens for cell ",
        tokenized$cell_names[[cell_index]], ".")
      if (pool_type == "all") {
        final_one <- hidden$select(2L, seq_len)$squeeze(1L)
        final_two <- hidden$select(2L, seq_len - 1L)$squeeze(1L)
        genes <- hidden$narrow(2L, 1L, seq_len - 2L)
        maximum <- torch::torch_max(genes, dim = 2L)[[1L]]$squeeze(1L)
        average <- torch::torch_mean(genes, dim = 2L)$squeeze(1L)
        pooled <- torch::torch_cat(list(final_one, final_two, maximum, average), dim = 1L)
      } else {
        pooled <- torch::torch_max(hidden, dim = 2L)[[1L]]$squeeze(1L)
      }
      embeddings[cell_index, ] <- as.numeric(torch::as_array(pooled$cpu()))
      rm(value_tensor, position_tensor, hidden, pooled)
    }
    if (isTRUE(verbose) && (batch_no %% progress_step == 0L || batch_no == length(starts))) {
      message(sprintf("scFoundation embedded %d/%d cells", max(index), n_cells))
    }
    if (batch_no %% 10L == 0L) invisible(gc(verbose = FALSE))
  }
  rownames(embeddings) <- tokenized$cell_names
  colnames(embeddings) <- paste0("scFoundation_", seq_len(ncol(embeddings)))
  object <- fmverse_to_seurat_v5(object)
  if (!is.null(cells)) object <- subset(object, cells = tokenized$cell_names)
  object[["scFoundation"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings, key = "scFoundation_", assay = tokenized$assay
  )
  structure(list(model_name = "scFoundation", embeddings = embeddings, seurat = object,
    model = model_module, audit = loaded$audit,
    inputs = list(values = lapply(tokenized$cells, `[[`, "values"),
      gene_positions = lapply(tokenized$cells, `[[`, "gene_positions"),
      cell_names = tokenized$cell_names, gene_order = tokenized$gene_order,
      assay = tokenized$assay, token_counts = tokenized$token_counts),
    analyses = list(), plots = list(), history = list(),
    metadata = list(checkpoint_name = resolved$name, checkpoint = resolved$path,
      config = config, device = device, batch_size = as.integer(batch_size), n_cells = n_cells,
      pool_type = pool_type, embedding_layer = "final xTrimoGene encoder output",
      pooling = if (pool_type == "all") "last total token + previous total token + gene max + gene mean" else "max across compacted gene and total tokens",
      matched_features = tokenized$matched_features, total_features = tokenized$total_features,
      token_counts = tokenized$token_counts, preprocessing = tokenized$preprocessing)),
    class = c("fmverse_result", "list"))
}
