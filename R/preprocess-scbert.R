# scBERT input adapter and R-native embedding inference.

fmverse_scbert_gene_index <- function(gene_index_file = NULL) {
  if (is.null(gene_index_file)) {
    roots <- c(file.path(getwd(), "00scBERT", "data"),
      file.path(getwd(), "..", "00scBERT", "data"))
    candidates <- file.path(roots, c("panglao_gene_index.tsv", "panglao_genes.tsv",
      "gene_index.tsv"))
    found <- candidates[file.exists(candidates)]
    if (!length(found)) stop("scBERT requires the official 16,906-gene Panglao order. ",
      "Extract Panglao `var_names` to 00scBERT/data/panglao_gene_index.tsv; see 00scBERT/README.md.",
      call. = FALSE)
    gene_index_file <- found[[1L]]
  }
  if (!file.exists(gene_index_file)) stop("scBERT gene index file does not exist: ",
    gene_index_file, call. = FALSE)
  ext <- tolower(tools::file_ext(gene_index_file))
  if (ext %in% c("csv")) {
    tab <- utils::read.csv(gene_index_file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    tab <- utils::read.delim(gene_index_file, stringsAsFactors = FALSE,
      check.names = FALSE, quote = "", comment.char = "")
  }
  gene <- if ("gene_name" %in% names(tab)) tab$gene_name else if ("gene" %in% names(tab)) {
    tab$gene
  } else tab[[1L]]
  gene <- trimws(as.character(gene))
  gene <- gene[!is.na(gene) & nzchar(gene)]
  if (length(gene) != 16906L || anyDuplicated(gene)) {
    stop("scBERT gene index must contain exactly 16,906 unique nonempty Panglao gene names.",
      call. = FALSE)
  }
  gene
}

fmverse_scbert_tokenize <- function(object, gene_index_file = NULL, assay = NULL,
                                    cells = NULL) {
  gene_order <- fmverse_scbert_gene_index(gene_index_file)
  object <- fmverse_to_seurat_v5(object)
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
  if (!assay %in% SeuratObject::Assays(object)) stop("Assay not found: ", assay)
  if (!"counts" %in% SeuratObject::Layers(object[[assay]])) {
    stop("scBERT input requires a raw counts layer.")
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
  gene_names <- if ("scGPT_gene_name" %in% colnames(feature_metadata)) {
    as.character(feature_metadata$scGPT_gene_name)
  } else rownames(counts)
  row_index <- match(gene_order, gene_names)
  matched <- which(!is.na(row_index))
  if (!length(matched)) stop("No Seurat features match the scBERT Panglao gene index.")
  matched_rows <- row_index[matched]
  # The official preprocessing keeps the first match for each reference symbol.
  keep <- !duplicated(matched_rows)
  matched <- matched[keep]
  matched_rows <- matched_rows[keep]
  aligned <- counts[matched_rows, , drop = FALSE]
  cell_total <- Matrix::colSums(aligned)
  if (any(!is.finite(cell_total) | cell_total <= 0)) {
    bad <- cell_names[!is.finite(cell_total) | cell_total <= 0]
    stop("scBERT normalization found cells with zero mapped counts: ",
      paste(head(bad, 5L), collapse = ", "))
  }
  # Scanpy normalize_total(target_sum=1e4), log1p(base=2), then the official
  # Dataset casts to long integers and clips the expression class at five.
  token_i <- integer(); token_j <- integer(); token_x <- integer()
  aligned <- methods::as(aligned, "dgCMatrix")
  for (j in seq_along(cell_names)) {
    first <- aligned@p[[j]] + 1L
    last <- aligned@p[[j + 1L]]
    if (first > last) next
    pos <- seq.int(first, last)
    counts_j <- aligned@x[pos]
    positive <- counts_j > 0
    if (!any(positive)) next
    # `aligned` contains only genes present in the query, so its sparse row
    # indices are compact. Restore each row's original Panglao position before
    # assigning the token; otherwise missing genes shift every later position.
    gene_pos <- matched[aligned@i[pos][positive] + 1L]
    log_expr <- log1p(counts_j[positive] / cell_total[[j]] * 1e4) / log(2)
    bins <- pmin(5L, pmax(0L, as.integer(log_expr)))
    nonzero_bin <- bins > 0L
    token_i <- c(token_i, rep.int(j, sum(nonzero_bin)))
    token_j <- c(token_j, gene_pos[nonzero_bin])
    token_x <- c(token_x, bins[nonzero_bin])
  }
  # The final sequence position is the official appended zero token.
  tokens <- Matrix::sparseMatrix(i = token_i, j = token_j, x = token_x,
    dims = c(length(cell_names), 16907L), giveCsparse = TRUE)
  rownames(tokens) <- cell_names
  list(gene_ids = tokens, cell_names = cell_names, gene_order = gene_order,
    assay = assay, matched_features = length(matched), total_features = nrow(counts),
    token_counts = rep.int(16907L, length(cell_names)),
    preprocessing = list(
      method = "scBERT official Panglao alignment; Scanpy normalize_total(1e4); log1p(base=2); integer truncation and upper clipping at 5",
      fixed_gene_order = normalizePath(gene_index_file, winslash = "/", mustWork = TRUE),
      gene_count = 16906L, unmatched_features_removed = nrow(counts) - length(matched),
      duplicate_input_symbols = sum(duplicated(gene_names)),
      appended_zero_token = TRUE
    ))
}

#' Embed cells with the R-native scBERT checkpoint.
#' @param object Seurat v4/v5 or `fmobj`.
#' @param model Model option, `scBERT`.
#' @param checkpoint `panglao-human` or local checkpoint file.
#' @param checkpoint_dir Optional checkpoint root.
#' @param gene_index_file Optional ordered Panglao gene-name file.
#' @param batch_size Cells per forward pass; keep small because input has 16,907 tokens.
#' @param cells Optional cell barcodes.
#' @param assay Counts assay.
#' @param device Torch device.
#' @param verbose Report inference progress.
#' @return An `fmverse_result` suitable for conversion to an `fmobj` run.
#' @export
.fmverse_scbert_legacy <- function(object, model = "scBERT", checkpoint = "panglao-human",
                                   checkpoint_dir = NULL, gene_index_file = NULL,
                                   batch_size = 1L, cells = NULL, assay = NULL,
                                   device = "cpu", verbose = interactive()) {
  if (length(batch_size) != 1L || is.na(batch_size) || batch_size < 1L ||
      batch_size != as.integer(batch_size)) stop("`batch_size` must be a positive integer.")
  resolved <- fmverse_resolve_scbert_checkpoint(checkpoint, checkpoint_dir)
  tokenized <- fmverse_scbert_tokenize(object, gene_index_file, assay, cells)
  model_module <- fmverse_scbert_build_model()
  loaded <- fmverse_load_scbert_weights(model_module, resolved$path)
  model_module <- loaded$model
  torch_device <- torch::torch_device(device)
  model_module$to(device = torch_device); model_module$eval()
  n_cells <- length(tokenized$cell_names)
  embeddings <- matrix(NA_real_, nrow = n_cells, ncol = 200L)
  starts <- seq.int(1L, n_cells, by = as.integer(batch_size))
  progress_step <- max(1L, ceiling(length(starts) / 20L))
  for (batch_no in seq_along(starts)) {
    first <- starts[[batch_no]]
    index <- first:min(first + as.integer(batch_size) - 1L, n_cells)
    input <- as.matrix(tokenized$gene_ids[index, , drop = FALSE])
    input_ids <- torch::torch_tensor(input, dtype = torch::torch_long(), device = torch_device)
    hidden <- torch::with_no_grad(model_module(input_ids, return_encodings = TRUE))
    pooled <- torch::torch_mean(hidden[, 1:16906, ], dim = 2L)
    embeddings[index, ] <- matrix(torch::as_array(pooled$cpu()), nrow = length(index), ncol = 200L)
    if (isTRUE(verbose) && (batch_no %% progress_step == 0L || batch_no == length(starts))) {
      message(sprintf("scBERT embedded %d/%d cells", max(index), n_cells))
    }
    rm(input_ids, hidden, pooled)
    if (batch_no %% 10L == 0L) invisible(gc(verbose = FALSE))
  }
  rownames(embeddings) <- tokenized$cell_names
  colnames(embeddings) <- paste0("scBERT_", seq_len(ncol(embeddings)))
  object <- fmverse_to_seurat_v5(object)
  if (!is.null(cells)) object <- subset(object, cells = tokenized$cell_names)
  object[["scBERT"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings, key = "scBERT_", assay = tokenized$assay
  )
  structure(list(model_name = "scBERT", embeddings = embeddings, seurat = object,
    model = model_module, audit = loaded$audit,
    inputs = list(gene_ids = tokenized$gene_ids, cell_names = tokenized$cell_names,
      gene_order = tokenized$gene_order, assay = tokenized$assay),
    analyses = list(), plots = list(), history = list(),
    metadata = list(checkpoint_name = resolved$name, checkpoint = resolved$path,
      device = device, batch_size = as.integer(batch_size), n_cells = n_cells,
      embedding_layer = "final Performer hidden states",
      pooling = "mean across the 16,906 Panglao gene positions; appended sentinel excluded",
      matched_features = tokenized$matched_features, total_features = tokenized$total_features,
      preprocessing = tokenized$preprocessing)),
    class = c("fmverse_result", "list"))
}
